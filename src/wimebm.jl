module WIMEBM # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import Integrals as Intgr
import NonlinearSolve as NlinSol

export Spectrum
export bretschneider

struct Spectrum
    freq::Vec
    period::Vec
    density::Vec

    Spectrum(freq::Vec, density::Vec) =
        length(freq) == length(density) ?
            new(freq, 2pi ./ freq, density) :
            throw(ArgumentError("Frequency and density vectors must be of the same length."))
end # struct Spectrum

function bretschneider(Hs::Float64, Tp::Float64, freq::Vec)::Spectrum
    T = 2pi ./ freq
    return Spectrum(freq, @. 1.25 * Hs^2 * T^5 / (8pi * Tp^4) * exp(-1.25(T/Tp)^4))
end # function bretschneider

function dispersion_relation(
    k::Complex{Float64}, p::@NamedTuple{omega::Float64, gamma::Float64, h::Float64, par::Par}
)::Complex{Float64}
    F = p.par.Y * p.h^3 / 12(1 - p.par.nu^2)
    lhs = (F*k^4 + p.par.rho * (p.par.g - 0.9p.h*p.omega^2) - im*p.omega*p.gamma) * k
    rhs = p.par.rho * p.omega^2
    return lhs - rhs
end # function dispersion_relation

function wavenumber_ice(omega::Float64, h::Float64, par::Par, gamma::Float64=0.0)::Complex{Float64}
    sol = NlinSol.solve(
        NlinSol.NonlinearProblem(
            dispersion_relation, omega^2/par.g + 0im, (; omega, gamma, h, par)
        ),
        NlinSol.NewtonRaphson(; autodiff=NlinSol.AutoFiniteDiff()),
    ) # solve
    if !NlinSol.SciMLBase.successful_retcode(sol)
        @warn "Nonlinear solver did not converge. Result may be inaccurate."
        @isdebugging() && @show (sol.retcode, sol.resid)
    end # if !
    return sol.u
end # function wavenumber_ice

function moment(S::Spectrum, n::Int; coeff::Function=one)::Float64
    int = Intgr.solve(
        Intgr.SampledIntegralProblem(@.(coeff(S.freq) * S.freq^n * S.density), S.freq),
        Intgr.SimpsonsRule()
    )
    if !Intgr.SciMLBase.successful_retcode(int)
        @warn "Integral did not converge. Result may be inaccurate."
        @isdebugging() && @show int.retcode
    end # if !
    return int.u
end # function moment

moment_elevation(S::Spectrum, n::Int)::Float64 = moment(S, n)

moment_strain(S::Spectrum, n::Int, h::Float64, par::Par)::Float64 = moment(
    S, n; coeff=(freq -> real(wavenumber_ice(freq, h, par, 0.0))^4 * h^2/4)
)

# TODO inplace version?
function attenuate(S::Spectrum, L::Float64, h::Float64, phi::Float64, par::Par)::Spectrum
    alpha = 2phi * imag.(wavenumber_ice.(S.freq, h, Ref(par), par.Gamma))
    return Spectrum(S.freq, @.(S.density * exp(-alpha*L)))
end # function attenuate

wave_period(S::Spectrum)::Float64 = 2pi * sqrt(moment_elevation(S, 0) / moment_elevation(S, 2))

function wave_length(S::Spectrum, h::Float64, par::Par)::Float64
    Tw = wave_period(S)
    kw = wavenumber_ice(2pi/Tw, h, par, 0.0)
    return 2pi / kw
end # function wave_length

wave_strain(S::Spectrum, h::Float64, par::Par)::Float64 = 2sqrt(moment_strain(S, 0, h, par))

fsd(d::Float64, dmn::Float64, dmx::Float64, par::Par)::Float64 =
    dmn <= d <= dmx ?
        par.gamma * dmn^par.gamma * d^(-(par.gamma+1)) / (1 - (dmn/dmx)^par.gamma) :
        0.0

function mean_size(dmn::Float64, dmx::Float64, par::Par)::Float64
    sol = Intgr.solve(
        Intgr.IntegralProblem((d, _) -> d * fsd(d, dmn, dmx, par), (dmn, dmx)),
        Intgr.QuadGKJL()
    )
    if !Intgr.SciMLBase.successful_retcode(sol)
        @warn "Integral did not converge. Result may be inaccurate."
        @isdebugging() && @show (sol.retcode, sol.resid)
    end # if !
    return sol.u
end # function mean_size


end # module WIMEBM
