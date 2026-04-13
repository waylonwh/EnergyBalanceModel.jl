module WIMEBM # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import Integrals as Intgr
import NonlinearSolve as NlinSol

export Spectrum
export bretschneider, monochromatic

struct Spectrum
    freq::Vec
    period::Vec
    density::Vec
end # struct Spectrum

Spectrum(freq::Vec, density::Vec) = length(freq) == length(density) ?
    Spectrum(freq, 2pi ./ freq, density) :
    throw(ArgumentError("Frequency and density vectors must be of the same length."))

Base.copy(S::Spectrum)::Spectrum = Spectrum(S.freq, S.period, S.density)

function bretschneider(
    Hs::Float64, Tp::Float64, freq::Vec=collect(range(2pi/23.8, 2pi/2.5; step=7.5e-2))
)::Spectrum
    T = 2pi ./ freq
    return Spectrum(freq, @. 1.25 * Hs^2 * T^5 / (8pi * Tp^4) * exp(-1.25(T/Tp)^4))
end # function bretschneider

monochromatic(A::Float64, T::Float64)::Spectrum = Spectrum(
    [2pi/T-1e-2, 2pi/T, 2pi/T+1e-2], [0.0, A^2/2, 0.0]
)

function dispersion_relation(
    k::ComplexF64, omega::Float64, gamma::Float64, h::Float64, par::Par
)::ComplexF64
    F = par.Y * h^3 / 12(1 - par.nu^2)
    lhs = (F*k^4 + par.rhow * (par.g - 0.9h*omega^2) - im*omega*gamma) * k
    rhs = par.rhow * omega^2
    return lhs - rhs
end # function dispersion_relation

function wavenumber_ice(
    omega::Float64, h::Float64, par::Par, gamma::Float64=0.0;
    init::ComplexF64=omega^2/par.g + 0im, abstol::Float64=1e-10
)::ComplexF64
    prob = NlinSol.NonlinearProblem(
        (k, _) -> dispersion_relation(k, omega, gamma, h, par), init
    )
    sol = NlinSol.solve(
        prob, NlinSol.NewtonRaphson(; autodiff=NlinSol.AutoFiniteDiff()); abstol
    )
    (real(sol.u) < 0 || imag(sol.u) < 0 || abs(sol.resid) > 1) && abs(init) < 100  &&
        return wavenumber_ice(omega, h, par, gamma; init=10init) # try another initial guess
    stalledsol = sol.retcode === NlinSol.ReturnCode.Stalled && abs(sol.resid) < 10abstol
    NlinSol.SciMLBase.successful_retcode(sol) || stalledsol ||
        @warn(
            "Nonlinear solver did not converge when solving for wavenumber. Result may be inaccurate.",
            sol.retcode, sol.resid
        )
    @isdebugging() && stalledsol && (
        isdefined(Main, :stalled_kice_sols) ?
            Main.stalled_kice_sols += 1 : @eval Main stalled_kice_sols = 1
    )
    return sol.u
end # function wavenumber_ice

function moment(S::Spectrum, n::Int; coeff::Function=one)::Float64
    int = Intgr.solve(
        Intgr.SampledIntegralProblem(@.(coeff(S.freq) * S.freq^n * S.density), S.freq),
        Intgr.SimpsonsRule()
    )
    if !Intgr.SciMLBase.successful_retcode(int)
        @warn "Integral did not converge when computing moment. Result may be inaccurate."
        @isdebugging() && @show int.retcode
    end # if !
    return int.u
end # function moment

moment_elevation(S::Spectrum, n::Int)::Float64 = moment(S, n)

moment_strain(S::Spectrum, n::Int, h::Float64, par::Par)::Float64 = moment(
    S, n; coeff=(freq -> real(wavenumber_ice(freq, h, par, 0.0))^4 * h^2/4)
)

@persistent(
    alpha1::Vector{Float64}, input::UInt=0x0,
    function attenuate(S::Spectrum, L::Float64, h::Float64, phi::Float64, par::Par)::Spectrum
        if input != hash((S.freq, h, par))
            alpha1 = imag.(wavenumber_ice.(S.freq, h, Ref(par), par.Gamma))
            input = hash((S.freq, h, par))
        end # if !=
        return Spectrum(S.freq, S.period, @. S.density * exp(-2phi*alpha1 * L))
    end # function attenuate
) # @persistent

wave_period(S::Spectrum)::Float64 = 2pi * sqrt(moment_elevation(S, 0) / moment_elevation(S, 2))

function wave_length(S::Spectrum, h::Float64, par::Par)::Float64
    Tw = wave_period(S)
    kw = wavenumber_ice(2pi/Tw, h, par, 0.0)
    return 2pi / kw
end # function wave_length

wave_strain(S::Spectrum, h::Float64, par::Par)::Float64 = 2sqrt(moment_strain(S, 0, h, par))

function fracture_distance(S::Spectrum, h::Float64, phi::Float64, L::Float64, par::Par)::Float64
    sol = NlinSol.solve(
        NlinSol.IntervalNonlinearProblem(
            (l, _) -> wave_strain(attenuate(S, l, h, phi, par), h, par), (0.0, L)
        )
    )
    if !NlinSol.SciMLBase.successful_retcode(sol)
        @warn "Nonlinear solver did not converge when solving for fracture distance. Result may be inaccurate." # TODO get rid of this warning
        @isdebugging() && nothing #@show (sol.retcode, sol.resid)
    end # if !
    return sol.u
end # function fracture_distance

fsd(d::Float64, dmn::Float64, dmx::Float64, par::Par)::Float64 =
    dmn <= d <= dmx ?
        par.gamma * dmn^par.gamma * d^(-(par.gamma+1)) / (1 - (dmn/dmx)^par.gamma) :
        0.0

function mean_size(dmx::Float64, par::Par)::Float64
    sol = Intgr.solve(
        Intgr.IntegralProblem((d, _) -> d * fsd(d, par.dmn, dmx, par), (par.dmn, dmx)),
        Intgr.QuadGKJL()
    )
    if !Intgr.SciMLBase.successful_retcode(sol)
        @warn "Integral did not converge when computing mean size. Result may be inaccurate."
        @isdebugging() && @show (sol.retcode, sol.resid)
    end # if !
    return sol.u
end # function mean_size

# Physical grid length at a given index in metres
function grid_length(st::SpaceTime{F}, i::Int)::Float64 where F
    if i == 1 # first grid
        dtheta = 1/2 * (asin(st.x[1]) + asin(st.x[2]))
    elseif i == st.nx # last grid
        dtheta = pi/2 - 1/2 * (asin(st.x[st.nx-1]) + asin(st.x[st.nx]))
    else # middle grids
        dtheta = 1/2 * (asin(st.x[i+1]) - asin(st.x[i-1]))
    end
    return dtheta / (pi/2) * 1e7
end # function grid_length

function Infrastructure.step!(
    ::WIModel, t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime{F}, par::Par; spectrum::Spectrum
)::Collection{Vec} where F
    Infrastructure.step!(MIZModel(), t, f, vars, st, par) # thermodynamics
    edgeinx = findfirst(>(0), vars.h)
    isnothing(edgeinx) && return vars # no ice
    xi = edgeinx # to retain scope of xi
    local S = copy(spectrum) # protect original spectrum
    local lastS = copy(S) # retain scope for partial breakup
    for outer xi in edgeinx:st.nx
        L = grid_length(st, xi)
        lastS = copy(S)
        S = attenuate(S, L, vars.h[xi], vars.phi[xi], par)
        wave_strain(S, vars.h[xi], par) < par.Ec && break # wave can not break ice anymore
        dmx = 1/2 * wave_length(S, vars.h[xi], par)
        dbar = mean_size(dmx, par)
        dbar < vars.D[xi] && (vars.D[xi] = dbar) # update floe size if it has been reduced by breaking # !
    end # for outer xi
    # partial breakup
    L = grid_length(st, xi)
    l = fracture_distance(lastS, vars.h[xi], vars.phi[xi], L, par)
    S = attenuate(lastS, l, vars.h[xi], vars.phi[xi], par)
    frontd = mean_size(1/2*wave_length(spectrum, vars.h[xi], par), par)
    dbar = (l*frontd + (L-l)*vars.D[xi]) / L
    dbar < vars.D[xi] && (vars.D[xi] = dbar) # !
    return vars
end # function Infrastructure.step!

end # module WIMEBM
