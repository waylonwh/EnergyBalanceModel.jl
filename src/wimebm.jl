# TODO docs & and docs for each submodule
module WIMEBM # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import EnergyBalanceModel.MIZEBM

import Integrals as Intgr
import NonlinearSolve as NlinSol

export Spectrum
export bretschneider, monochromatic

struct Spectrum{T<:AbstractFloat} <: AbstractSpectrum
    freq::Vector{T}
    period::Vector{T}
    density::Vector{T}
end # struct Spectrum

function Spectrum(freq::Vector{T}, density::Vector{T}) where {T<:AbstractFloat}
    return length(freq) == length(density) ?
        Spectrum{T}(freq, T(2pi) ./ freq, density) :
        throw(ArgumentError("Frequency and density vectors must be of the same length."))
end

function bretschneider(Hs::T, Tp::T, freq::Vector{T})::Spectrum{T} where {T<:AbstractFloat}
    period = T(2pi) ./ freq
    return Spectrum(freq, @. T(1.25) * Hs^2 * period^5 / (T(8pi) * Tp^4) * exp(-T(1.25) * (period/Tp)^4))
end # function bretschneider

function bretschneider(Hs::T, Tp::T)::Spectrum{T} where {T<:AbstractFloat}
    freq = collect(range(T(2pi)/T(23.8), T(2pi)/T(2.5); step=T(7.5e-2)))
    return bretschneider(Hs, Tp, freq)
end

function monochromatic(Hs::T, Tp::T, freq::Vector{T}; eps::Real=1e-6) where {T<:AbstractFloat}
    epsT = T(eps)
    return Spectrum(freq, @. Hs^2 / T(16) * exp(-(freq - T(2pi)/Tp)^2 / (T(2)*epsT)) / sqrt(T(2pi) * epsT))
end

function monochromatic(Hs::T, Tp::T; eps::Real=1e-6) where {T<:AbstractFloat}
    epsT = T(eps)
    freq = collect(range(T(2pi)/(Tp+T(0.1)), T(2pi)/(Tp-T(0.1)); step=T(1e-3)))
    return monochromatic(Hs, Tp, freq; eps=epsT)
end

function dispersion_relation(
    k::Complex{T}, omega::T, gamma::T, h::T, par::Par{T}
)::Complex{T} where {T<:AbstractFloat}
    F = par.Y * h^3 / (T(12) * (one(T) - par.nu^2))
    lhs = (F*k^4 + par.rhow * (par.g - T(0.9)*h*omega^2) - im*omega*gamma) * k
    rhs = par.rhow * omega^2
    return lhs - rhs
end # function dispersion_relation

function wavenumber_ice(
    omega::T, h::T, par::Par{T}, gamma::T=zero(T);
    init::Union{Nothing,Complex{T}}=nothing, abstol::Real=1e-10
)::Complex{T} where {T<:AbstractFloat}
    _init = isnothing(init) ? complex(omega^2/par.g, zero(T)) : init
    _abstol = T(abstol)
    prob = NlinSol.NonlinearProblem{false}(
        (k, _) -> dispersion_relation(k, omega, gamma, h, par), _init
    )
    sol = NlinSol.solve(
        prob, NlinSol.NewtonRaphson(; autodiff=NlinSol.AutoFiniteDiff()); abstol=_abstol
    )
    (real(sol.u) < zero(T) || imag(sol.u) < zero(T) || abs(sol.resid) > one(T)) && abs(_init) < T(100) &&
        return wavenumber_ice(omega, h, par, gamma; init=T(10)*_init, abstol=_abstol) # try another initial guess
    stalledsol = sol.retcode === NlinSol.ReturnCode.Stalled && abs(sol.resid) < T(10)*_abstol
    NlinSol.SciMLBase.successful_retcode(sol) || stalledsol ||
        @warn(
            "Nonlinear solver did not converge when solving for wavenumber. Result may be inaccurate.",
            sol.retcode, sol.resid
        )
    @isdebugging() && stalledsol && (
        isdefined(@__MODULE__, :_stalled_kice_sols) ?
            (global _stalled_kice_sols += 1) : (global _stalled_kice_sols = 1)
    )
    return sol.u
end # function wavenumber_ice

function moment(S::Spectrum{T}, n::Int; coeff::Function=one)::T where {T<:AbstractFloat}
    int = Intgr.solve(
        Intgr.SampledIntegralProblem(@.(coeff(S.freq) * S.freq^n * S.density), S.freq),
        Intgr.SimpsonsRule()
    )
    Intgr.SciMLBase.successful_retcode(int) ||
        @warn "Integral did not converge when computing moment. Result may be inaccurate." int.retcode
    return int.u
end # function moment

function moment_elevation(S::Spectrum{T}, n::Int)::T where {T<:AbstractFloat}
    return moment(S, n)
end

function moment_strain(S::Spectrum{T}, n::Int, h::T, par::Par{T})::T where {T<:AbstractFloat}
    return moment(S, n; coeff=(freq -> real(wavenumber_ice(freq, h, par, zero(T)))^4 * h^2/T(4)))
end

@persistent(
    alpha1 = Float64[],
    input::UInt=0x0,
    function attenuate(S::Spectrum, L::Real, h::Real, phi::Real, par::Par)::Spectrum
        T = eltype(S.freq)
        if input != hash((S.freq, h, par))
            alpha1 = imag.(wavenumber_ice.(S.freq, h, Ref(par), par.Gamma))
            input = hash((S.freq, h, par))
        end # if !=
        return Spectrum{T}(S.freq, S.period, @. S.density * exp(-T(2)*phi*alpha1 * L))
    end # function attenuate
) # @persistent

function wave_period(S::Spectrum{T})::T where {T<:AbstractFloat}
    return T(2pi) * sqrt(moment_elevation(S, 0) / moment_elevation(S, 2))
end

function wave_length(S::Spectrum{T}, h::T, par::Par{T})::T where {T<:AbstractFloat}
    Tw = wave_period(S)
    kw = wavenumber_ice(T(2pi)/Tw, h, par, zero(T))
    return real(T(2pi) / kw)
end # function wave_length

function wave_strain(S::Spectrum{T}, h::T, par::Par{T})::T where {T<:AbstractFloat}
    return T(2) * sqrt(moment_strain(S, 0, h, par))
end

function fracture_distance(S::Spectrum{T}, h::T, phi::T, L::T, par::Par{T})::T where {T<:AbstractFloat}
    sol = NlinSol.solve(
        NlinSol.IntervalNonlinearProblem(
            (l, _) -> wave_strain(attenuate(S, l, h, phi, par), h, par) - par.Ec, (zero(T), L)
        )
    )
    NlinSol.SciMLBase.successful_retcode(sol) ||
        @warn(
            "Nonlinear solver did not converge when solving for fracture distance. Result may be inaccurate.",
            sol.retcode, sol.resid
        )
    return sol.u
end # function fracture_distance

function fsd(d::T, dmn::T, dmx::T, par::Par{T})::T where {T<:AbstractFloat}
    return dmn <= d <= dmx ?
        par.gamma * dmn^par.gamma * d^(-(par.gamma+1)) / (1 - (dmn/dmx)^par.gamma) :
        zero(T)
end

function mean_size(dmx::T, par::Par{T})::T where {T<:AbstractFloat}
    sol = Intgr.solve(
        Intgr.IntegralProblem((d, _) -> d * fsd(d, par.dmn, dmx, par), (par.dmn, dmx)),
        Intgr.QuadGKJL()
    )
    Intgr.SciMLBase.successful_retcode(sol) ||
        @warn(
            "Integral did not converge when computing mean size. Result may be inaccurate.",
            sol.retcode, sol.resid
        )
    return sol.u
end # function mean_size

# Physical grid length at a given index in metres
function grid_length(st::SpaceTime{F,T}, i::Int)::T where {F,T<:AbstractFloat}
    if i == 1 # first grid
        dtheta = one(T)/2 * (asin(st.x[1]) + asin(st.x[2]))
    elseif i == st.nx # last grid
        dtheta = T(pi)*(1//2) - one(T)*(1//2) * (asin(st.x[st.nx-1]) + asin(st.x[st.nx]))
    else # middle grids
        dtheta = one(T)/2 * (asin(st.x[i+1]) - asin(st.x[i-1]))
    end
    return dtheta / (T(pi)*(1//2)) * T(1e7)
end # function grid_length

function updateD!(
    newD::T, xi::Int, vars::Collection{Vector{T}}, l::Real=1, L::Real=1
)::Nothing where {T<:AbstractFloat}
    _l = T(l)
    _L = T(L)
    dbar = (_l*newD + (_L-_l)*vars.D[xi]) / _L # weighted average based on fracture distance
    dbar < vars.D[xi] && (vars.D[xi] = dbar) # update floe size if it has been reduced by breaking # !
    return nothing
end # function updateD!

function Infrastructure.initialise(
    model::WIModel, st::SpaceTime{F,T}, forcing::Forcing{T,C}, par::Par{T}, init::Collection{Vector{T}};
    lastonly::Bool=true, spectrum::Spectrum{T}
) where {F,C,T<:AbstractFloat} # -> Tuple{Collection{Vector}, Solutions{WIModel,F,V}, Solutions{WIModel,F,V}}
    vars, sols, annusol = MIZEBM._initialise(model, st, forcing, par, init; lastonly)
    sols.spectrum_ref[] = deepcopy(spectrum) # store spectrum in sols for later reference
    annusol.spectrum_ref[] = deepcopy(spectrum)
    vars.Ewave = zeros(T, st.nx)
    vars.lambda = Vector{T}(undef, st.nx)
    return (vars, sols, annusol)
end # function Infrastructure.initialise

function Infrastructure.step!(
    ::WIModel,
    t::T,
    f::T,
    vars::Collection{Vector{T}},
    st::SpaceTime{F,T},
    par::Par{T};
    spectrum::Spectrum{T}
)::Collection{Vector{T}} where {F,T<:AbstractFloat}
    Infrastructure.step!(MIZModel(), t, f, vars, st, par) # thermodynamics
    edgeinx = findfirst(>(zero(T)), vars.h)
    if isnothing(edgeinx) # no ice
        vars.Ewave .= zero(T)
        vars.lambda .= wave_period(spectrum)
        return vars
    elseif edgeinx > 1 # at least once cell has no ice
        vars.Ewave[1:edgeinx-1] .= zero(T)
        vars.lambda[1:edgeinx-1] .= wave_period(spectrum)
    end # if isnothing, elseif
    # attenuate spectrum
    spect = spectrum
    for xi in edgeinx:st.nx
        L = grid_length(st, xi)
        atted_spect = attenuate(spect, L, vars.h[xi], vars.phi[xi], par)
        atted_strain = wave_strain(atted_spect, vars.h[xi], par)
        half_atted_spect = attenuate(atted_spect, L/T(2), vars.h[xi], vars.phi[xi], par)
        half_atted_lambda = wave_length(half_atted_spect, vars.h[xi], par)
        if atted_strain > par.Ec # full breakup
            dbar = mean_size(one(T)/2*half_atted_lambda, par)
            updateD!(dbar, xi, vars)
        elseif wave_strain(spect, vars.h[xi], par) > par.Ec # partial breakup
            l = fracture_distance(spect, vars.h[xi], vars.phi[xi], L, par)
            half_partial_atted_spect = attenuate(spect, l/T(2), vars.h[xi], vars.phi[xi], par)
            half_partial_atted_lambda = wave_length(half_partial_atted_spect, vars.h[xi], par)
            frontd = mean_size(one(T)/2*half_partial_atted_lambda, par)
            updateD!(frontd, xi, vars, l, L)
        end # if >, else
        spect = atted_spect # update spectrum
        vars.Ewave[xi] = wave_strain(half_atted_spect, vars.h[xi], par)
        vars.lambda[xi] = half_atted_lambda
    end # for xi
    return vars
end # function Infrastructure.step!

end # module WIMEBM
