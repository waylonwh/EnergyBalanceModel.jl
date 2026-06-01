module WIMEBM # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import EnergyBalanceModel.MIZEBM

import InteractiveUtils as IU
import Integrals as Intgr
import NonlinearSolve as NlinSol

export Spectrum
export bretschneider, monochromatic

"""
    Spectrum(freq::Vec, density::Vec)

Represents a wave energy spectrum for the wave-ice interaction model ([`WIModel`](@ref)).

The angular frequencies `freq` (rad s⁻¹) and the corresponding spectral energy `density`
must be vectors of the same length. The wave periods are computed as `2π ./ freq` and stored
in the `period` field.

# Fields
- `freq::Vec`: angular frequencies (rad s⁻¹)
- `period::Vec`: wave periods (s), derived as `2π ./ freq`
- `density::Vec`: spectral energy density at each frequency

See also [`bretschneider`](@ref) and [`monochromatic`](@ref) for convenience constructors.

# Examples
```julia-repl
julia> Spectrum([0.5, 1.0, 1.5], [0.2, 0.5, 0.1])
Spectrum([0.5, 1.0, 1.5], [12.566370614359172, 6.283185307179586, 4.1887902047863905], [0.2, 0.5, 0.1])
```
"""
struct Spectrum <: AbstractSpectrum
    freq::Vec
    period::Vec
    density::Vec
end # struct Spectrum

Spectrum(freq::Vec, density::Vec) = length(freq) == length(density) ?
    Spectrum(freq, 2pi ./ freq, density) :
    throw(ArgumentError("Frequency and density vectors must be of the same length."))

Base.show(io::IO, S::Spectrum)::Nothing = print(
    io, "Spectrum(", length(S.freq), " components)"
)

function Base.show(io::IO, ::MIME"text/plain", S::Spectrum)::Nothing
    println(io, "Spectrum with ", length(S.freq), " frequency components:")
    println(io, "  ω ∈ [", round(minimum(S.freq); digits=3), ", ", round(maximum(S.freq); digits=3), "] rad s⁻¹")
    println(io, "  T ∈ [", round(minimum(S.period); digits=2), ", ", round(maximum(S.period); digits=2), "] s")
    print(io, "  peak density at T = ", round(S.period[argmax(S.density)]; digits=2), " s")
    return nothing
end # function Base.show

struct WavenumberCache{T<:AbstractFloat}
    key::UInt # hash of (T, freqs, Gamma, abstal, par)
    dh::T
    hmax::T
    wavenumber::Matrix{Complex{T}} # freqs × hs
end

const _wavenumber_ice_cache_ref = Ref{NTuple{2,WavenumberCache}}()

"""
    bretschneider(Hs::Float64, Tp::Float64, freq::Vec=collect(range(2π/23.8, 2π/2.5; step=7.5e-2))) -> Spectrum

Construct a Bretschneider wave [`Spectrum`](@ref) with significant wave height `Hs` (m) and
peak period `Tp` (s), evaluated at the angular frequencies `freq` (rad s⁻¹).

The spectral energy density is given by
``S(T) = \\frac{1.25 H_s^2 T^5}{8\\pi T_p^4} \\exp\\!\\left[-1.25 (T/T_p)^4\\right]``,
where ``T = 2\\pi / \\omega`` is the wave period.

# Examples
```julia-repl
julia> S = bretschneider(3.0, 9.5);

julia> length(S.freq)
14
```
"""
function bretschneider(
    Hs::Float64, Tp::Float64, freq::Vec=collect(range(2pi/23.8, 2pi/2.5; step=7.5e-2))
)::Spectrum
    T = 2pi ./ freq
    return Spectrum(freq, @. 1.25 * Hs^2 * T^5 / (8pi * Tp^4) * exp(-1.25(T/Tp)^4))
end # function bretschneider

"""
    monochromatic(Hs::Float64, Tp::Float64, freq::Vec=collect(range(2π/(Tp+0.1), 2π/(Tp-0.1); step=1e-3)); eps::Float64=1e-6) -> Spectrum

Construct an approximately monochromatic wave [`Spectrum`](@ref) with significant wave
height `Hs` (m) and peak period `Tp` (s).

The energy is concentrated near the peak angular frequency ``2\\pi / T_p`` using a narrow
Gaussian of variance `eps`, so that the total energy matches that of a monochromatic wave of
height `Hs`. Smaller `eps` produces a sharper peak.

# Examples
```julia-repl
julia> S = monochromatic(3.0, 9.5);

julia> isapprox(S.freq[argmax(S.density)], 2pi / 9.5; atol=1e-3)
true
```
"""
monochromatic(
    Hs::Float64, Tp::Float64, freq::Vec=collect(range(2pi/(Tp+0.1), 2pi/(Tp-0.1); step=1e-3));
    eps::Float64=1e-6
)::Spectrum = Spectrum(freq, @. Hs^2 / 16 * exp(-(freq - 2pi/Tp)^2 / 2eps) / sqrt(2pi * eps))

function dispersion_relation(
    k::ComplexF64, omega::Float64, gamma::Float64, h::Float64, par::Par
)::ComplexF64
    F = par.Y * h^3 / 12(1 - par.nu^2)
    lhs = (F*k^4 + par.rhow * (par.g - 0.9h*omega^2) - im*omega*gamma) * k
    rhs = par.rhow * omega^2
    return lhs - rhs
end # function dispersion_relation

function get_cache(
    freqs::AbstractVector, h::T, gamma::T, abstol::T, par::Collection{T}
)::WavenumberCache{T} where T <: AbstractFloat
    index = iszero(gamma) ? 1 : 2
    cache = _wavenumber_ice_cache_ref[][index]::WavenumberCache{T}
    if hash((T, freqs, gamma, abstol, par)) != cache.key || h > cache.hmax + 1
        @warn "Wavenumber cache miss. Recomputing."
        cache = cache_wavenumber!(freqs, par, 1//10^4, max(10, h+1); abstol)[index]
    end # if ||
    return cache
end # function get_cache

function interpolate_wavenumber(h::AbstractFloat, cache::WavenumberCache) # -> Union{Nothing,Vector{Complex}}
    h > cache.hmax && ArgumentError("h is out of bounds for wavenumber cache.")
    col = floor(Int, h / cache.dh) + 1
    weight = h % cache.dh / cache.dh
    interp = weight * cache.wavenumber[:,col]
    @. interp += (1-weight) * cache.wavenumber[:,col+1]
    return interp
end # function interpolate_wavenumber

function wavenumber_ice(
    omega::Float64, h::Float64, par::Par, gamma::Float64=0.0;
    abstol::Float64=1e-10, # init::ComplexF64=omega^2/par.g + 0im
)::ComplexF64
    prob = NlinSol.NonlinearProblem{false}(
        (k, _) -> dispersion_relation(k, omega, gamma, h, par), omega^2/par.g + 0im
    )
    sol = NlinSol.solve(
        prob, NlinSol.NewtonRaphson(; autodiff=NlinSol.AutoFiniteDiff()); abstol
    )
    # (real(sol.u) < 0 || imag(sol.u) < 0 || abs(sol.resid) > 1) && abs(init) < 100  &&
    #     return wavenumber_ice(omega, h, par, gamma; init=10init) # try another initial guess
    stalledsol = sol.retcode === NlinSol.ReturnCode.Stalled && abs(sol.resid) < 10abstol
    NlinSol.SciMLBase.successful_retcode(sol) || # stalledsol ||
        @warn(
            "Nonlinear solver did not converge when solving for wavenumber. Result may be inaccurate.",
            sol.retcode, sol.resid, sol.u, omega, h
        )
    @isdebugging() && stalledsol && (
        isdefined(@__MODULE__, :_stalled_kice_sols) ?
            (global _stalled_kice_sols += 1) : (global _stalled_kice_sols = 1)
    )
    return sol.u
end # function wavenumber_ice

function wavenumber_ice(
    omegas::AbstractVector, h::T, par::Collection{T}, gamma::T=zero(T); abstol::T=T(1e-10)
)::Vector{Complex{T}} where T <: AbstractFloat
    cache = get_cache(omegas, h, gamma, abstol, par)
    return h < cache.hmax ?
        interpolate_wavenumber(h, cache) :
        wavenumber_ice.(omegas, h, Ref(par), gamma; abstol)
end # function wavenumber_ice

function cache_wavenumber!(
    freqs::AbstractVector, par::Collection{T}, dh::Real, hmax::Real; abstol::T=T(1e-10)
)::NTuple{2,WavenumberCache{T}} where T <: AbstractFloat
    hvec = Vector{T}(0:dh:hmax)
    tup = (
        WavenumberCache(
            hash((T, freqs, zero(T), abstol, par)), T(dh), T(hmax),
            wavenumber_ice.(freqs, hvec', Ref(par), zero(T); abstol)
        ),
        WavenumberCache(
            hash((T, freqs, par.Gamma, abstol, par)), T(dh), T(hmax),
            wavenumber_ice.(freqs, hvec', Ref(par), par.Gamma; abstol))
    )
    _wavenumber_ice_cache_ref[] = tup
    return tup
end # function cache_wavenumber

function moment(S::Spectrum, n::Int; coeff::Vector=ones(length(S.freq)))::Float64
    int = Intgr.solve(
        Intgr.SampledIntegralProblem(@.(coeff * S.freq^n * S.density), S.freq),
        Intgr.SimpsonsRule()
    )
    Intgr.SciMLBase.successful_retcode(int) ||
        @warn "Integral did not converge when computing moment. Result may be inaccurate." int.retcode
    return int.u
end # function moment

moment_elevation(S::Spectrum, n::Int)::Float64 = moment(S, n)

moment_strain(S::Spectrum, n::Int, h::Float64, par::Par)::Float64 = moment(
    S, n; coeff=real.(wavenumber_ice(S.freq, h, par, 0.0)).^4 * h^2/4
)

function attenuate(S::Spectrum, L::Float64, h::Float64, phi::Float64, par::Par)::Spectrum
    alpha1 = imag.(wavenumber_ice(S.freq, h, par, par.Gamma))
    return Spectrum(S.freq, S.period, @. S.density * exp(-2phi*alpha1 * L))
end # function attenuate) # @persistent

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
            (l, _) -> wave_strain(attenuate(S, l, h, phi, par), h, par) - par.Ec, (0.0, L)
        )
    )
    NlinSol.SciMLBase.successful_retcode(sol) ||
        @warn(
            "Nonlinear solver did not converge when solving for fracture distance. Result may be inaccurate.",
            sol.retcode, sol.resid
        )
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
    Intgr.SciMLBase.successful_retcode(sol) ||
        @warn(
            "Integral did not converge when computing mean size. Result may be inaccurate.",
            sol.retcode, sol.resid
        )
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

function updateD!(newD::Float64, xi::Int, vars::Collection{Vec}, l::Float64=1.0, L::Float64=1.0)::Nothing
    dbar = (l*newD + (L-l)*vars.D[xi]) / L # weighted average based on fracture distance
    dbar < vars.D[xi] && (vars.D[xi] = dbar) # update floe size if it has been reduced by breaking # !
    return nothing
end # function updateD!

function Infrastructure.initialise(
    model::WIModel, st::SpaceTime, forcing::Forcing, par::Par, init::Collection{Vec};
    lastonly::Bool=true, spectrum::Spectrum
) # -> Tuple{Collection{Vec}, Solutions{WIModel,F,V}, Solutions{WIModel,F,V}}
    vars, sols, annusol = MIZEBM._initialise(model, st, forcing, par, init; lastonly)
    sols.spectrum_ref[] = deepcopy(spectrum) # store spectrum in sols for later reference
    annusol.spectrum_ref[] = deepcopy(spectrum)
    vars.Ewave = zeros(st.nx)
    vars.lambda = Vec(undef, st.nx)
    cache_wavenumber!(spectrum.freq, par, 1//10^4, 10)
    return (vars, sols, annusol)
end # function Infrastructure.initialise

function Infrastructure.step!(
    ::WIModel, t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime, par::Par; spectrum::Spectrum
)::Collection{Vec}
    breakup = falses(st.nx) # track which cells are breaking
    edgeinx = findfirst(>(0), vars.h)
    if isnothing(edgeinx) # no ice
        vars.Ewave .= 0.0
        vars.lambda .= wave_length(spectrum, 0.0, par)
        edgeinx = st.nx + 1 # skip loop
    elseif edgeinx > 1 # at least once cell has no ice
        vars.Ewave[1:edgeinx-1] .= 0.0
        vars.lambda[1:edgeinx-1] .= wave_length(spectrum, 0.0, par)
    end # if isnothing, elseif
    # attenuate spectrum
    spect = spectrum
    for xi in edgeinx:st.nx
        L = grid_length(st, xi)
        atted_spect = attenuate(spect, L, vars.h[xi], vars.phi[xi], par)
        atted_strain = wave_strain(atted_spect, vars.h[xi], par)
        half_atted_spect = attenuate(spect, L/2, vars.h[xi], vars.phi[xi], par)
        half_atted_lambda = wave_length(half_atted_spect, vars.h[xi], par)
        if atted_strain > par.Ec # full breakup
            dbar = mean_size(1/2*half_atted_lambda, par)
            updateD!(dbar, xi, vars)
            breakup[xi] = true
        elseif wave_strain(spect, vars.h[xi], par) > par.Ec # partial breakup
            l = fracture_distance(spect, vars.h[xi], vars.phi[xi], L, par)
            half_partial_atted_spect = attenuate(spect, l/2, vars.h[xi], vars.phi[xi], par)
            half_partial_atted_lambda = wave_length(half_partial_atted_spect, vars.h[xi], par)
            frontd = mean_size(1/2*half_partial_atted_lambda, par)
            updateD!(frontd, xi, vars, l, L)
            breakup[xi] = true
        end # if >, else
        spect = atted_spect # update spectrum
        vars.Ewave[xi] = wave_strain(half_atted_spect, vars.h[xi], par)
        vars.lambda[xi] = half_atted_lambda
    end # for xi
    Infrastructure.step!(MIZModel(), t, f, vars, st, par; breakup) # thermodynamics
    return vars
end # function Infrastructure.step!

end # module WIMEBM
