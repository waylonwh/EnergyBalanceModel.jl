module Infrastructure # EnergyBalanceModel.

using ..Utilities

import Integrals as Intgr, InteractiveUtils as IU, SparseArrays as SA, Statistics as Stats, StyledStrings as SS

export AbstractModel, ClassicModel, MIZModel, WIModel, ModelDiff
export Collection, Forcing, Par, Solutions, SpaceTime, Vec, EBMProblem
export Spectrum, bretschneider, monochromatic
export default_parameters, default_parval
export get_diffop
export hemispheric_mean, ice_area
export create_storages, integrate, solve

"""
    Vec

Alias for `Vector{Float64}` to represent model variables.
"""
const Vec = Vector{Float64}

"""
    AbstractModel

Abstract type for the different energy balance models.
"""
abstract type AbstractModel end

"""
    MIZModel <: AbstractModel

Singleton type representing the extended idealised climate model with a marginal ice zone
(MIZ).
"""
struct MIZModel <: AbstractModel end

"""
    WIModel <: AbstractModel

Singleton type representing the wave ice interaction model, as an extension of the
[`MIZModel`](@ref).
"""
struct WIModel <: AbstractModel end

"""
    ClassicModel <: AbstractModel

Singleton type representing the classic idealised climate model by Wagner & Eisenman (2015).
"""
struct ClassicModel <: AbstractModel end

"""
    ModelDiff{B<:AbstractModel, S<:AbstractModel} <: AbstractModel

A type representing the difference (S-B) between two models `B` and `S`.
"""
struct ModelDiff{B<:AbstractModel, S<:AbstractModel} <: AbstractModel end

"""
    Collection{V}(args...)

A simple wrapper around `Dict{Symbol,V}` to allow dot syntax access to fields. Use syntax
for constructing a `Dict{Symbol,V}` to create a `Collection{V}`.

# Examples
```julia-repl
julia> parameters = Collection{Float64}(:D => 0.6, :A => 193.0, :B => 2.1)
Collection{Float64} with 3 entries:
  :A => 193.0
  :D => 0.6
  :B => 2.1

julia> parameters.D
0.6

julia> getproperty(parameters, :A)
193.0

julia> parameters[:A]
193.0

julia> parameters.F = 0.0; parameters.F
0.0
```
"""
struct Collection{V} <: AbstractDict{Symbol,V}
    dict::Dict{Symbol,V}
    Collection{V}(args...) where V = new{V}(Dict{Symbol,V}(args...))
end # struct Collection

Base.getproperty(coll::Collection, key::Symbol) = getindex(getfield(coll, :dict), key) # -> V
Base.setproperty!(coll::Collection, key::Symbol, val) = setindex!(getfield(coll, :dict), val, key) # -> Dict{Symbol,V}
Base.propertynames(coll::Collection)::Set{Symbol} = Set(keys(getfield(coll, :dict)))
Base.getindex(coll::Collection, key::Symbol) = getproperty(coll, key) # -> V
Base.setindex!(coll::Collection, val, key::Symbol) = setproperty!(coll, key, val) # -> Dict{Symbol,V}
Base.keys(coll::Collection)::Set{Symbol} = propertynames(coll)
Base.iterate(coll::Collection) = iterate(getfield(coll, :dict)) # -> Tuple{Pair{Symbol,V},Int} or Nothing
Base.iterate(coll::Collection, state::Int) = iterate(getfield(coll, :dict), state) # -> Tuple{Pair{Symbol,V},Int} or Nothing
Base.length(coll::Collection)::Int = length(getfield(coll, :dict))
Base.hash(coll::Collection, h::UInt)::UInt = hash(getfield(coll, :dict), h)

function uniqueunion(ca::Collection{A}, cb::Collection{B}) where {A, B}
    vtype = typejoin(A, B)
    vtype === Any && (vtype = Union{A, B})
    overlap = intersect(propertynames(ca), propertynames(cb))
    return all(key -> ca[key] === cb[key], overlap) ?
        Collection{vtype}() : Collection{vtype}(union(ca, cb))
end # function uniqueunion

"""
    Par

Alias for `Collection{Float64}` to represent model parameters.

See also [`Collection`](@ref).
"""
const Par = Collection{Float64}

"""
    SpaceTime{F}(urange::NTuple{2,Float64}, nx::Int, nt::Int, dur::Int; winter::Float64=0.26125, summer::Float64=0.77375)

Defines the spatial and temporal grid for the model.

Type parameter `F` is a function that maps the uniform grid `u` with `nx` gridboxes in the
range `urange` to the model grid `x`. `urange` and `F` should be chosen such that
F(urange)=[0,1] and F(u₂)>F(u₁) for u₂>u₁.

`nt` is the number of timesteps per year, and `dur` is the duration of the simulation in
years. `winter` and `summer` are the times in a year (in [0,1]) when winter and summer peaks
occur.

    SpaceTime(nx::Int, nt::Int, dur::Int; kwargs...)
    SpaceTime{identity}(nx::Int, nt::Int, dur::Int; kwargs...)

A convenience constructor for `SpaceTime{identity}` with `urange=(0.0, 1.0)`.

    SpaceTime{sin}(nx::Int, nt::Int, dur::Int; kwargs...)

A convenience constructor for `SpaceTime{sin}` with `urange=(0.0, π/2)`.

# Examples
```julia-repl
julia> SpaceTime(100, 2000, 30)
SpaceTime{identity} with:
  100 latitudinal gridboxes: [0.005, 0.015, 0.025, 0 … 5, 0.975, 0.985, 0.995]
  2000 timesteps per year: [0.00025, 0.00075, 0.001 … 99875, 0.99925, 0.99975]
  30 years of simulation: t∈[0,30]
  winter at t=0.26125, summer at t=0.77375

julia> SpaceTime{sin}(180, 2000, 30)
SpaceTime{sin} with:
  180 latitudinal gridboxes: [0.00436331, 0.0130896, … 762, 0.999914, 0.99999]
  2000 timesteps per year: [0.00025, 0.00075, 0.001 … 99875, 0.99925, 0.99975]
  30 years of simulation: t∈[0,30]
  winter at t=0.26125, summer at t=0.77375
```
"""
struct SpaceTime{F}
    nx::Int # number of evenly spaced latitudinal gridboxes (equator to pole)
    u::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64} # grid before scale
    x::Vec # grid
    dur::Int # duration of simulation in years
    nt::Int # number of timesteps per year (limited by numerical stability)
    dt::Float64 # timestep
    t::Vec # time vector in a year
    T::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64} # full time series
    winter::@NamedTuple{t::Float64, inx::Int}
    summer::@NamedTuple{t::Float64, inx::Int}

    function SpaceTime{F}(
        urange::NTuple{2,Float64}, nx::Int, nt::Int, dur::Int;
        winter::Float64=0.26125, summer::Float64=0.77375
    ) where F
        du = (urange[2]-urange[1]) / nx
        u = range(urange[1] + du/2, urange[2] - du/2, nx)
        x = F.(u)
        dt = 1 / nt
        t = collect(range(dt/2, 1 - dt/2, nt))
        T = dt/2 : dt : dur - dt/2
        winterinx = round(Int, nt*winter)
        summerinx = round(Int, nt*summer)
        return new{F}(
            nx, u, x, dur, nt, dt, t, T, (t=winter, inx=winterinx), (t=summer, inx=summerinx)
        )
    end # function SpaceTime{F}
end # struct SpaceTime{F}

SpaceTime{identity}(nx::Int, nt::Int, dur::Int; kwargs...) = SpaceTime{identity}((0.0, 1.0), nx, nt, dur; kwargs...)
SpaceTime{sin}(nx::Int, nt::Int, dur::Int; kwargs...) = SpaceTime{sin}((0.0, pi/2), nx, nt, dur; kwargs...)
SpaceTime(args...; kwargs...) = SpaceTime{identity}(args...; kwargs...)

Base.show(io::IO, st::SpaceTime)::Nothing = print(
    io, typeof(st), '(', st.nx, ", ", st.nt, ", ", st.dur, ')'
)

function Base.show(io::IO, ::MIME"text/plain", st::SpaceTime)::Nothing
    println(io, typeof(st), " with:")

    nxstr = "  $(st.nx) latitudinal gridboxes: "
    buffer = iobuffer(io)
    show(buffer, st.x)
    vecstr = ctruncate(String(take!(buffer.io)), displaysize(io)[2]-length(nxstr)-2, " … ")
    println(io, nxstr, vecstr)

    nystr = "  $(st.nt) timesteps per year: "
    buffer = iobuffer(io)
    show(buffer, st.t)
    vecstr = ctruncate(String(take!(buffer.io)), displaysize(io)[2]-length(nystr)-2, " … ")
    println(io, nystr, vecstr)

    println(io, "  $(st.dur) years of simulation: t∈[0,$(st.dur)]")
    print(io, "  winter at t=$(st.winter.t), summer at t=$(st.summer.t)")
    return nothing
end # function Base.show

"""
    Forcing(base::Float64)

Defines a constant climate forcing of value `base`.

    Forcing(base::Float64, peak::Float64, cool::Float64, holdyrs::NTuple{2,Int}, rates::NTuple{2,Float64})

Defines a time-varying climate forcing that first holds at `base` for `holdyrs[1]` years,
warms to `peak` at rate `rates[1]>0` per year, holds at `peak` for `holdyrs[2]` years, cools
to `cool` at rate `rates[2]<0` per year, and then holds at `cool` thereafter. Warming and
cooling times (`(peak-base)/rates[1]` and `(cool-peak)/rates[2]`) must be positive integers.

`Forcing` has a field `domain` which is a tuple of years at which the forcing pattern
changes: `Forcing.domain=(1:0, 2:start of warming, 3:reach peak, 4:start of cooling,
5:reach cool)`.

A `Forcing` object can be called as a function to evaluate the forcing at a given time in
years.

# Examples
```julia-repl
julia> Forcing(0.0)
Forcing{true}(0.0) is constant:
  F(t)=0.0, t∈[0,∞)

julia> f = Forcing(0.0, 5.0, -5.0, (10, 10), (0.5, -0.5))
Forcing{false} varies from 0.0 up to 5.0 and back to -5.0:
  F(t)={  0.0            , t∈[ 0,10) (base)
       {  0.0 + 0.5(t-10), t∈[10,20) (warming)
       {  5.0            , t∈[20,30) (peak)
       {  5.0 - 0.5(t-30), t∈[30,50) (cooling)
       { -5.0            , t∈[50, ∞) (cool)

julia> f.domain
(0, 10, 20, 30, 50)

julia> f(17.57)
3.785
```
"""
struct Forcing{V}
    base::Float64 # base forcing
    peak::Float64 # peak forcing
    cool::Float64 # forcing after cooldown
    holdyrs::NTuple{2,Int} # years to hold at (base, peak) forcing
    rates::NTuple{2,Float64} # rates of change
    domain::NTuple{5,Int} # years at which forcing pattern changes

    # constant forcing
    Forcing(base::Float64) = new{false}(
        base, base, base, (0, 0), (0.0, 0.0), (0, 0, 0, 0, 0)
    )
    # warming/cooling forcing
    function Forcing(
        base::Float64, peak::Float64, cool::Float64, holdyrs::NTuple{2,Int}, rates::NTuple{2,Float64}
    )
        domainvec = zeros(Int, 5)
        # hold at base
        @. domainvec[2:5] += holdyrs[1]
        # time to warm
        warming = (peak-base) / rates[1]
        rates[1]>0 && isinteger(warming) ?
            @.(domainvec[3:5] += warming) :
            throw(ArgumentError("Warming time must be positive integer. Got $warming y."))
        # hold at peak
        @. domainvec[4:5] += holdyrs[2]
        # time to cool
        cooling = (cool-peak) / rates[2]
        rates[2]<0 && isinteger(cooling) ?
            domainvec[5] += cooling :
            throw(ArgumentError("Cooling time must be positive integer. Got $cooling y."))
        return new{true}(base, peak, cool, holdyrs, rates, Tuple(domainvec))
    end # function Forcing
end # struct Forcing{F}

function Base.show(io::IO, forcing::Forcing{false})::Nothing
    print(io, "Forcing(", forcing.base, ')')
    printstyled(io, " (constant forcing)", color=:light_black)
    return nothing
end # function Base.show

Base.show(io::IO, forcing::Forcing{true})::Nothing = print(
    io,
    "Forcing(", forcing.base, ", ", forcing.peak, ", ", forcing.cool, ')'
)

function Base.show(io::IO, ::MIME"text/plain", forcing::Forcing{false})::Nothing
    println(io, typeof(forcing), '(', forcing.base, ") is constant:")
    print(io, "  F(t)=", forcing.base, ", t∈[0,∞)")
    return nothing
end # function Base.show

function Base.show(io::IO, ::MIME"text/plain", forcing::Forcing{true})::Nothing
    println(
        io,
        typeof(forcing), " varies from ", forcing.base, " up to ", forcing.peak, " and back to ", forcing.cool, ':'
    )
    head = "  F(t)={ "
    headpad = lpad("{ ", length(head))
    biaslen = maximum(length∘string, (forcing.base, forcing.peak, forcing.cool))
    ratelen = maximum(length∘string∘abs, forcing.rates)
    domainlen = maximum(length∘string, forcing.domain)
    constline(field::Symbol)::String = string(
        lpad(getfield(forcing, field), biaslen), " "^(ratelen+domainlen+7)
    )
    varyline(bias::Float64, rate::Float64, start::Int)::String = string(
        lpad(bias, biaslen),
        ' ', rate>0 ? '+' : '-', ' ',
        lpad(abs(rate), ratelen), "(t-", lpad(start, domainlen), ")"
    )
    domainstr(i::Int, nextdomain::String=string(forcing.domain[i+1]))::String = string(
        ", t∈[", lpad(forcing.domain[i], domainlen), ',', lpad(nextdomain, domainlen), ")"
    )
    print(io, head, constline(:base), domainstr(1))
    printstyled(io, " (base)\n", color=:light_black)
    print(io, headpad, varyline(forcing.base, forcing.rates[1], forcing.domain[2]), domainstr(2))
    printstyled(io, " (warming)\n", color=:light_black)
    print(io, headpad, constline(:peak), domainstr(3))
    printstyled(io, " (peak)\n", color=:light_black)
    print(io, headpad, varyline(forcing.peak, forcing.rates[2], forcing.domain[4]), domainstr(4))
    printstyled(io, " (cooling)\n", color=:light_black)
    print(io, headpad, constline(:cool), domainstr(5, "∞"))
    printstyled(io, " (cool)", color=:light_black)
end # function Base.show

# evaluate forcing at time T (in years)
(forcing::Forcing{false})(::Float64)::Float64 = forcing.base # constant forcing
function (forcing::Forcing{true})(T::Float64)::Float64 # varying forcing
    if T < forcing.domain[2] # hold at base
        return forcing.base
    elseif T < forcing.domain[3] # warming
        return forcing.base + forcing.rates[1] * (T-forcing.domain[2])
    elseif T < forcing.domain[4] # hold at peak
        return forcing.peak
    elseif T < forcing.domain[5] # cooling
        return forcing.peak + forcing.rates[2] * (T-forcing.domain[4])
    else # hold at cool
        return forcing.cool
    end # if <, elseif*3, else
end # function (forcing::Forcing{true})

# Spectrum for WIM
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
struct Spectrum
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

"""
    Solutions{M,F,V}

An object to store model solutions. Type parameter `M` is the model type (`MIZModel` or
`ClassicModel`); `F` is the function used to map the uniform grid to the model grid in
`SpaceTime{F}`; `V` is a boolean indicating whether the climate forcing is variable.
`V` is `true` for variable forcing.

# Fields
- `model::M`: model type
- `spacetime::SpaceTime{F}`: space and time on which solutions are defined
- `ts::Vec`: time vector for stored solutions
- `forcing::Forcing{V}`: climate forcing
- `parameters::Par`: model parameters
- `initconds::Collection{Vec}`: initial conditions
- `lastonly::Bool`: whether to store solutions for each time step only for the last year
- `raw::Collection{Vector{Vec}}`: solutions for each time step
- `annual::@NamedTuple{winter::..., summer::..., avg::...}`: seasonal peak and annual
    average solution storage for each year

Note that each value of `raw`, `annual.winter`, `annual.summer`, and `annual.avg` is
a vector of vectors. For example, `raw.E[ti]::Vector{Float64}` stores the solution for
enthalpy at time step `ts[ti]::Float64`, and `annual.avg.T[y]::Vector{Float64}` stores
the annual average temperature for year `y::Int`.
"""
struct Solutions{M<:AbstractModel,F,V}
    spacetime::SpaceTime{F} # space and time which solutions are defined on
    ts::Vec # time vector for stored solution
    forcing::Forcing{V} # climate forcing
    parameters::Par # model parameters
    initconds::Collection{Vec} # initial conditions
    lastonly::Bool # store only last year of solution
    raw::Collection{Vector{Vec}} # solution storage
    annual::@NamedTuple{
        winter::Collection{Vector{Vec}}, summer::Collection{Vector{Vec}}, avg::Collection{Vector{Vec}}
    } # seasonal peak and annual avg
    spectrum_ref::Ref{Spectrum} # spectrum used for WIModel, if applicable

    function Solutions{M}(
        st::SpaceTime{F}, forcing::Forcing{V}, par::Par, init::Collection{Vec},
        vars::Set{Symbol}, lastonly::Bool=true
    ) where {M<:AbstractModel, F, V} # Solutions
        if lastonly
            dur_store = 1
            ts::Vec = st.dur-1 + st.dt/2 : st.dt : st.dur - st.dt/2
        else # !lastonly
            dur_store = st.dur
            ts = st.dt/2 : st.dt : st.dur - st.dt/2
        end # if lastonly, else
        # construct raw solution storage
        solraw = Collection{Vector{Vec}}()
        foreach(var -> (solraw[var] = Vector{Vec}(undef, length(ts))), vars)
        # construct seasonal solution storage template
        seasonaltemp = Collection{Vector{Vec}}()
        foreach(var -> (seasonaltemp[var] = Vector{Vec}(undef, st.dur)), vars)
        return new{M,F,V}(
            st, # spacetime
            ts,
            forcing,
            par, # parameters
            init, # initconds
            lastonly,
            solraw, # raw
            (
                winter=deepcopy(seasonaltemp),
                summer=deepcopy(seasonaltemp),
                avg=deepcopy(seasonaltemp)
            ), # ( # seasonal
            Ref{Spectrum}() # spectrum_ref
        ) # new
    end # function Solutions
end # struct Solutions{M,F,V}

function Base.:-(
    sx::Solutions{X,F,false}, sy::Solutions{Y,F,false}
)::Solutions{ModelDiff{X,Y},F,false} where {X<:AbstractModel, Y<:AbstractModel, F}
    (sx.spacetime.x == sy.spacetime.x && sx.spacetime.t == sy.spacetime.t) ||
        throw(
            ArgumentError(
                "Cannot compute difference of solutions defined on different space-time grids."
            )
        ) # throw
    st = sx.spacetime.dur == sy.spacetime.dur ?
        sx.spacetime :
        SpaceTime{F}(sx.spacetime.x, sx.spacetime.t, max(sx.spacetime.dur, sy.spacetime.dur))
    forcing = Forcing(sx.forcing.base - sy.forcing.base)
    par = uniqueunion(sx.parameters, sy.parameters)
    init = uniqueunion(sx.initconds, sy.initconds)
    vars = intersect(propertynames(sx.raw), propertynames(sy.raw))
    lastonly = sx.lastonly || sy.lastonly
    diffsol = Solutions{ModelDiff{X,Y}}(st, forcing, par, init, vars, lastonly)
    xinx = findall(in(diffsol.ts), sx.ts)
    yinx = findall(in(diffsol.ts), sy.ts)
    foreach(var -> (diffsol.raw[var] = sx.raw[var][xinx] .- sy.raw[var][yinx]), vars)
    for season in 1:3, var in vars
        diffsol.annual[season][var] = sx.annual[season][var][1:st.dur] .- sy.annual[season][var][1:st.dur]
    end # for season, var
    return diffsol
end # function Base.:-

Base.show(io::IO, sols::Solutions)::Nothing = print(
    io,
    typeof(sols), '(',
    sols.spacetime.nx, '×', length(sols.ts), "@(", first(sols.ts), ':', sols.spacetime.dt, ':', last(sols.ts), "), ",
    propertynames(sols.raw),
    ')'
)

function Base.show(io::IO, ::MIME"text/plain", sols::Solutions)::Nothing
    println(io, typeof(sols), " with:")
    println(io, "  ", length(sols.raw), " solution variables: ", propertynames(sols.raw))
    xhead = "  on $(sols.spacetime.nx) latitudinal gridboxes: "
    buffer = iobuffer(io)
    show(buffer, sols.spacetime.x)
    vecstr = ctruncate(String(take!(buffer.io)), displaysize(io)[2]-length(xhead)-2, " … ")
    println(io, xhead, vecstr)
    println(io, "  and " , length(sols.ts), " timesteps: ", first(sols.ts), ':', sols.spacetime.dt, ':', last(sols.ts))
    print(io, "  with forcing ", repr(sols.forcing))
    return nothing
end # function Base.show

get_spectrum(sol::Solutions{WIModel})::Spectrum = sol.spectrum_ref[] # -> Spectrum

const _secyear = 31536000 # number of seconds in a year

# default parameter values
const default_parval = Par(
    # classic
    :D => 0.6, # diffusivity for heat transport (W m^-2 K^-1)
    :A => 193.0, # OLR when T = T_m (W m^-2)
    :B => 2.1, # OLR temperature dependence (W m^-2 K^-1)
    :cw => 9.8, # ocean mixed layer heat capacity (W y m^-2 K^-1)
    :S0 => 420.0, # insolation at equator  (W m^-2)
    :S1 => 338.0, # insolation seasonal dependence (W m^-2)
    :S2 => 240.0, # insolation spatial dependence (W m^-2)
    :a0 => 0.7, # ice-free co-albedo at equator
    :a2 => 0.1, # ice-free co-albedo spatial dependence
    :ai => 0.4, # co-albedo where there is sea ice
    :Fb => 4.0, # heat flux from ocean below (W m^-2)
    :k => 2.0, # sea ice thermal conductivity (W m^-2 K^-1)
    :Lf => 9.5, # sea ice latent heat of fusion (W y m^-3)
    :cg => 1e-3 * 9.8, # ghost layer heat capacity(W y m^-2 K^-1)
    :tau => 1e-5 * 9.8, # ghost layer coupling timescale (y)
    # MIZ
    :Tm => 0.0, # mean temperature (C)
    :m1 => 1.6e-6_secyear, # empirical constants of lateral melt (m y^-1 K^-1)
    :m2 => 1.36, # empirical constants of lateral melt
    :alpha => 0.66, # floe geometry constant, Ai = alpha * D^2
    :rl => 0.5, # lead region width (m)
    :Dmin => 1.0, # new pancake size (m)
    :Dmax => 500.0, # largest floe length (m)
    :hmin => 0.1, # new pancake thickness (m)
    :kappa => 0.01_secyear, # floe welding parameter (m^-2 y^-1)
    :rhow => 1025.0, # Water density (kg m^-3)
    :cp => 3980.0 / _secyear, # Specific heat for seawater near freezing (W y kg^-1 K^-1)
    :ch => 6e-3, # Heat transfer coefficient
    :u0 => 0.01_secyear, # Surface friction velocity (m y^-1)
    # WIM
    :Y => 5.5e9, # Effective Young's modulus (Pa)
    :nu => 0.3, # Poisson's ratio
    :g => 9.81, # Gravitational acceleration (m s^-2),
    :Ec => 7.05e-5, # Breaking significant strain
    :Gamma => 13.0, # Viscous damping parameter (Pa m s^-1)
    :gamma => 2 + log2(0.9), # Power law exponent for floe size distribution
) # Par

# parameters used in each model
const classicmodel_parvars = Set{Symbol}(
    (:D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :cg, :tau)
)
const mizmodel_parvars = push!(
    copy(classicmodel_parvars),
    :Tm, :m1, :m2, :alpha, :rl, :Dmin, :Dmax, :hmin, :kappa, :rhow, :cp, :ch, :u0
)
const wimodel_parvars = push!(
    copy(mizmodel_parvars),
    :Y, :nu, :g, :Ec, :Gamma, :gamma
)

# Create a parameter dictionary from default values for a given Set
function default_parameters(paramset::Set{Symbol})::Par
    setvec = collect(paramset)
    return Par(setvec .=> getproperty.(Ref(default_parval), setvec))
end # function get_defaultparameters

"""
    default_parameters(<:AbstractModel) -> Par

Get default parameters for a given model.

# Examples
```julia-repl
julia> default_parameters(ClassicModel())
Collection{Float64} with 16 entries:
  :a2 => 0.1
  :F  => 0.0
  :A  => 193.0
  :k  => 2.0
  :D  => 0.6
  :S1 => 338.0
  :B  => 2.1
  :cw => 9.8
  :S2 => 240.0
  :S0 => 420.0
  ⋮   => ⋮
```
"""
function default_parameters end # stub
for model in IU.subtypes(AbstractModel)
    namelower = lowercase(split(string(model), '.')[end])
    @eval default_parameters(::$model)::Par = default_parameters($(Symbol(namelower, "_parvars")))
end # for model

(
    initconds_fromT(::ClassicModel, T::Vector{FT}, cw::FT,)::Collection{Vector{FT}}
) where FT<:AbstractFloat = Collection{Vector{FT}}(:Tg => T, :E => cw .* T)

(
    initconds_fromT(::AbstractModel, T::Vector{FT}, cw::FT)::Collection{Vector{FT}}
) where FT<:AbstractFloat = Collection{Vector{FT}}(
    :Tg => T, :Ei => zeros(FT, length(T)), :Ew => cw .* T,
    :h => zeros(FT, length(T)), :D => zeros(FT, length(T))
)

initvarset(::ClassicModel)::Set{Symbol} = Set((:Tg, :E))
initvarset(::AbstractModel)::Set{Symbol} = Set((:Tg, :Ei, :Ew, :h, :D))

function check_initconds(model::AbstractModel, initconds::Collection{<:Vector})::Nothing
    varset = initvarset(model)
    list(itr) = join(string.(itr), ", ", " and ")
    issubset(varset, propertynames(initconds)) ||
        throw(ArgumentError("Initial conditions for $(typeof(model)) must include $(list(varset))."))
    length(varset) == length(propertynames(initconds)) ||
        @warn "$(typeof(model)) only needs initial conditions for $(list(varset)). Extra variables will be ignored."
    return nothing
end # function check_initconds

"""
    EBMProblem(model::AbstractModel, ::Type{FT}=Float64; st, forcing, parameters, initconds, spectrum) -> EBMProblem{FT}

Bundle a `model` together with everything needed to integrate it — the space-time domain,
the climate forcing, the model parameters, the initial conditions, and (for `WIModel`) the
incident wave spectrum — into a single `EBMProblem` that can be passed to `solve`.

The floating-point type `FT` (default `Float64`) sets the precision used throughout the
problem. Every other input is supplied as a keyword argument; each one is optional and falls
back to a sensible default, so `EBMProblem(model)` produces a fully-specified, runnable
problem. The keyword arguments are:
- `st`: the `SpaceTime` domain, or a `NamedTuple` of overrides for its fields `nx`, `nt`,
  `dur`, and `F` (e.g. `st=(; nx=90, dur=100)`). Defaults to `SpaceTime{sin}(180, 2000,
  50)`.
- `forcing`: a `Forcing`, a `Number` (used as a constant forcing), or a `NamedTuple` of the
  fields accepted by `Forcing`. Defaults to no forcing, `Forcing(0)`.
- `parameters`: a `Collection` of parameters, or a `NamedTuple` of overrides applied on top
  of `default_parameters(model)`. Defaults to `default_parameters(model)`.
- `initconds`: a `Collection` of initial conditions, or a `Vector` giving a starting
  temperature field `T` (from which the remaining variables are derived). Defaults to a
  uniform temperature of 17°C with no ice.
- `spectrum`: a `Spectrum` describing the incident wave field. Only valid for `WIModel`,
  where it defaults to `bretschneider(3, 9.5)`; supplying it for any other model is an
  error.

# Examples
```julia-repl
julia> EBMProblem(WIModel())
EBMProblem{Float64} with:
  model:      WIModel
  spacetime:  SpaceTime{sin}(180, 2000, 50)
  forcing:    Forcing(0.0) (constant forcing)
  parameters: 32 entries
  initconds:  5 variables: Set([:Tg, :Ei, :D, :h, :Ew])
  spectrum:   Spectrum(30 components)

julia> EBMProblem(WIModel(); st=(dur=100,))
EBMProblem{Float64} with:
  model:      WIModel
  spacetime:  SpaceTime{sin}(180, 2000, 100)
  forcing:    Forcing(0.0) (constant forcing)
  parameters: 32 entries
  initconds:  5 variables: Set([:Tg, :Ei, :D, :h, :Ew])
  spectrum:   Spectrum(30 components)
```
"""
mutable struct EBMProblem{T<:AbstractFloat}
    model::AbstractModel
    st::SpaceTime
    forcing::Forcing
    parameters::Collection{T}
    initconds::Collection{Vec}
    spectrum::Union{Spectrum,Nothing}

    function EBMProblem{FT}(
        model::M,
        st::SpaceTime=SpaceTime{sin}(180, 2000, 50),
        forcing::Forcing=Forcing(zero(FT)),
        parameters::Collection{FT}=default_parameters(model),
        initconds::Union{Collection{Vector{FT}},Nothing}=nothing;
        spectrum::Union{Spectrum,Nothing}=nothing
    ) where {FT<:AbstractFloat, M<:AbstractModel}
        (M === AbstractModel || model isa ModelDiff) &&
            throw(
                ArgumentError("model must be one of the following types: ClassicModel, MIZModel, or WIModel.")
            )
        isnothing(initconds) ?
            initconds = initconds_fromT(model, fill(FT(17), st.nx), parameters.cw) :
            check_initconds(model, initconds)
        model isa WIModel || isnothing(spectrum) ||
            @warn "spectrum will be ignored for non-WIModel."
        model isa WIModel && isnothing(spectrum) &&
            (spectrum = bretschneider(FT(3), FT(9.5)))
        return new{FT}(model, st, forcing, parameters, initconds, spectrum)
    end # function EBMProblem{FT}
end # mutable struct EBMProblem

function EBMProblem(
    model::AbstractModel, ::Type{FT}=Float64;
    st::Union{SpaceTime,NamedTuple}=(; ),
    forcing::Union{Forcing,Number,NamedTuple}=0,
    parameters::Union{Collection{FT},NamedTuple}=(; ),
    initconds::Union{Collection{Vector{FT}},Vector{FT},Nothing}=nothing,
    spectrum::Union{Spectrum,Nothing}=nothing
)::EBMProblem{FT} where FT<:AbstractFloat
    if !(st isa SpaceTime)
        nx = get(st, :nx, 180)
        st_inst = SpaceTime{get(st, :F, sin)}(nx, get(st, :nt, 2000), get(st, :dur, 50))
    else # st isa SpaceTime
        nx = st.nx
        st_inst = st
    end # if !, else
    if forcing isa Forcing
        forcing_inst = forcing
    else # forcing isa Number or NamedTuple
        forcing_inst = forcing isa Number ?
            Forcing(FT(forcing)) :
            Forcing(forcing.base, forcing.peak, forcing.cool, forcing.holdyrs, forcing.rates)
    end # if isa, else
    if parameters isa Collection
        parameters_inst = parameters
    else # parameters isa NamedTuple
        parameters_inst = default_parameters(model)
        foreach(k -> (parameters_inst[k] = parameters[k]), propertynames(parameters))
    end # if isa, else
    if initconds isa Collection
        initconds_inst = initconds
    elseif initconds isa Vector
        initconds_inst = initconds_fromT(model, initconds, parameters_inst.cw)
    else # initconds isa nothing
        initconds_inst = initconds_fromT(model, fill(FT(17), nx), parameters_inst.cw)
    end # if isa, elseif, else
    return EBMProblem{FT}(model, st_inst, forcing_inst, parameters_inst, initconds_inst; spectrum)
end # function EBMProblem

Base.show(io::IO, prob::EBMProblem)::Nothing = print(
    io, typeof(prob), '(', typeof(prob.model), ", ", prob.st, ", ", repr(prob.forcing), ')'
)

function Base.show(io::IO, ::MIME"text/plain", prob::EBMProblem)::Nothing
    println(io, typeof(prob), " with:")
    println(io, "  model:      ", typeof(prob.model))
    println(io, "  spacetime:  ", prob.st)
    println(io, "  forcing:    ", repr(prob.forcing))
    println(io, "  parameters: ", length(prob.parameters), " entries")
    println(io, "  initconds:  ", length(prob.initconds), " variables: ", propertynames(prob.initconds))
    if !isnothing(prob.spectrum)
        print(io, "  spectrum:   ", prob.spectrum)
    else
        print(io, "  spectrum:   nothing")
    end # if !isnothing
    return nothing
end # function Base.show

# calculate diffusion operator matrix
@persistent(
    diffop::SA.SparseMatrixCSC{Float64,Int64} = SA.spzeros(Float64, 0, 0),
    @inline function get_diffop(st::SpaceTime{identity})::SA.SparseMatrixCSC{Float64,Int64}
        if size(diffop) != (st.nx, st.nx) # recalculate diffusion operator
            dx = 1 / st.nx
            xb = dx:dx:1-dx
            lambda = @. (1 - xb^2) / dx^2
            l1 = pushfirst!(-copy(lambda), 0.0)
            l2 = push!(-copy(lambda), 0.0)
            l3 = -l1 - l2
            diffop = SA.spdiagm(-1 => -l1[2:st.nx], 0 => -l3, 1 => -l2[1:st.nx-1])
        end # if !=
        return diffop
    end # function get_diffop
) # @persistent

@persistent(
    diffop::SA.SparseMatrixCSC{Float64,Int64} = SA.spzeros(Float64, 0, 0),
    xid::UInt = UInt(0),
    @inline function get_diffop(st::SpaceTime)::SA.SparseMatrixCSC{Float64,Int64}
        if xid != objectid(st.x) # recalculate diffusion operator
            x = [-st.x[1]; st.x; 2 - st.x[end]] # include ghost points
            diffx = diff(x)
            f = diffx[2:st.nx+1] # (xⱼ₊₁ - xⱼ) or Δⱼ
            b = -diffx[1:st.nx] # (xⱼ₋₁ - xⱼ) or ∇ⱼ
            l = 2:st.nx # lower diagonal row indices
            u = 1:st.nx-1 # upper diagonal row indices
            A1 = @. f / (b * (f-b))
            B1 = @. -(b+f) / (b*f)
            C1 = @. b / (f * (b-f))
            first = -2st.x .* SA.spdiagm(
                -1 => A1[l],
                0 => B1 + [A1[1]; zeros(Float64, st.nx-2); C1[st.nx]],
                1 => C1[u]
            ) # -2x ∂/∂x
            A2 = @. 2 / (b * (b-f))
            B2 = @. 2 / (b*f)
            C2 = @. 2 / (f * (f-b))
            second = (1 .- st.x.^2.0) .* SA.spdiagm( # `2.0` for v1.10 LTS compatibility
                -1 => A2[l],
                0 => B2 + [A2[1]; zeros(Float64, st.nx-2); C2[st.nx]],
                1 => C2[u]
            ) # (1-x^2) ∂²/∂x²
            diffop = first + second
            xid = objectid(st.x)
        end # if !=
        return diffop
    end # function get_diffop
) # @persistent

function annual_mean(annusol::Solutions)::Collection{Vec}
    # calculate annual mean for each variable except temperatures
    means = Collection{Vec}()
    for var in propertynames(annusol.raw)
        vecvec = annusol.raw[var]
        @boundscheck length(vecvec) == annusol.spacetime.nt ||
            throw(
                ArgumentError(
                    "Length of raw solution vector for $var does not match the number of timesteps per year, when calculating annual mean."
                )
            ) # throw
        means[var] = crossmean(vecvec)
    end # for var
    return means
end # function annual_mean

"""
    annual_mean(forcing::Forcing, st::SpaceTime, year::Int) -> Float64

Calculate the annual mean of the climate forcing for a given year.

# Examples
```julia-repl
julia> forcing = Forcing(0.0, 10.0, -5.0, (20, 10), (0.5, -0.5));

julia> st = SpaceTime{sin}(180, 2000, 30);

julia> annual_mean(forcing, st, 24)
1.75
```
"""
annual_mean(forcing::Forcing, st::SpaceTime, year::Int)::Float64 = Stats.mean(forcing.(year-1 .+ st.t))

function savesol!(
    sols::Solutions{M,F,C}, annusol::Solutions{M,F,C}, vars::Collection{Vec}, tinx::Int
)::Solutions{M,F,C} where {M<:AbstractModel, F, C}
    varscp = deepcopy(vars) # avoid reference issues
    year = ceil(Int, sols.spacetime.T[tinx])
    ti = mod1(tinx, sols.spacetime.nt) # index of time in the year
    # save raw data to annual
    foreach(var -> annusol.raw[var][ti] = varscp[var], propertynames(annusol.raw))
    # save raw data
    if !sols.lastonly # save all raw data
        foreach(var -> (sols.raw[var][tinx] = varscp[var]), propertynames(sols.raw))
    elseif tinx > length(sols.spacetime.T) - sols.spacetime.nt # save the raw data of the last year
        foreach(var -> (sols.raw[var][ti] = varscp[var]), propertynames(sols.raw))
    end # if !, elseif
    # save seasonal data
    if ti == sols.spacetime.winter.inx
        foreach(var -> (sols.annual.winter[var][year] = varscp[var]), propertynames(sols.annual.winter))
    elseif ti == sols.spacetime.summer.inx
        foreach(var -> (sols.annual.summer[var][year] = varscp[var]), propertynames(sols.annual.summer))
    elseif ti == sols.spacetime.nt # calculate annual average
        means = annual_mean(annusol)
        foreach(var -> (sols.annual.avg[var][year] = means[var]), propertynames(sols.annual.avg))
    end # if ==, elseif*2
    return sols
end # function savesol!

"""
    hemispheric_mean(vec::Vec, x::Vec) -> Float64

Calculate the hemispheric mean value of `vec` defined on grid `x` using the trapezoidal
rule.

# Examples
```julia-repl
julia> x = sin.(range(0, pi/2, 180));

julia> vec = @. 7.5 + 20(1 - 2x^2);

julia> hemispheric_mean(vec, x)
14.166324413879554
```
"""
function hemispheric_mean(vec::Vec, x::Vec)::Float64
    int = Intgr.solve(
        Intgr.SampledIntegralProblem(@.(2vec / (pi * sqrt(1-x^2))), x), Intgr.SimpsonsRule()
    )
    if !Intgr.SciMLBase.successful_retcode(int)
        @warn "Integral did not converge when computing hemispheric mean. Result may be inaccurate."
        @isdebugging() && @show int.retcode
    end # if !
    return int.u
end # function hemispheric_mean

"""
    ice_area(phi::Vec, x::Vec) -> Float64

Calculate the area covered by sea ice from the sea ice concentration `phi` defined on grid
`x` by discretised integration using the Simpson's rule.

# Examples
```julia-repl
julia> x = sin.(range(0, pi/2, 181))[1:end-1]; # avoid the singularity at x=1

julia> phi = @. 2x - x^2;

julia> ice_area(phi, x)
1.4075808945373096
```
"""
function ice_area(phi::Vec, x::Vec)::Float64
    int = Intgr.solve(
        Intgr.SampledIntegralProblem(@.(pi*x * phi / 2sqrt(1-x^2)), x), Intgr.SimpsonsRule()
    )
    if !Intgr.SciMLBase.successful_retcode(int)
        @warn "Integral did not converge when computing ice area. Result may be inaccurate."
        @isdebugging() && @show int.retcode
    end # if !
    return int.u
end # function ice_area

"""
    ice_area(sols::Solutions, season::Symbol, year::Int) -> Float64

Calculate the area covered by sea ice for a given season (must be `:winter`, `:summer`, or
`:avg`) and year from the solutions.

# Examples
```julia-repl
julia> sols = run_example();
[output omitted]

julia> ice_area(sols, :summer, 30)
0.43981792357403693
```
"""
ice_area(sols::Solutions{ClassicModel}, season::Symbol, year::Int)::Float64 =
    ice_area((getproperty(sols.annual, season).E[year].<0), sols.spacetime.x)
ice_area(sols::Solutions{MIZModel}, season::Symbol, year::Int)::Float64 =
    ice_area(getproperty(sols.annual, season).phi[year], sols.spacetime.x)

# stub for functions for each model
function step! end
function initialise end

function create_storages(
    ::M, solvars::Set{Symbol}, st::SpaceTime, forcing::Forcing, par::Par, init::Collection{Vec};
    lastonly::Bool
) where M<:AbstractModel # -> Tuple{Collection{Vec},Solutions{M,F,C},Solutions{M,F,C}}
    vars = deepcopy(init)
    sols = Solutions{M}(st, forcing, par, init, solvars, lastonly)
    annusol = Solutions{M}(st, forcing, par, init, solvars, true) # for calculating annual means
    return (vars, sols, annusol)
end # function create_storages

"""
    integrate(model::Union{MIZModel,ClassicModel}, st::SpaceTime, forcing::Forcing, par::Par, init::Collection{Vec}; lastonly::Bool=true, updatefreq::Float64=1.0) -> Solutions{ClassicModel,F,C}
    integrate(model::WIModel, st::SpaceTime, forcing::Forcing, par::Par, init::Collection{Vec}; lastonly::Bool=true, updatefreq::Float64=1.0, spectrum::Spectrum) -> Solutions{M,F,C}

Integrate the specified model over the given `SpaceTime` with climate `Forcing`, model
parameters `par`, and initial conditions `init`. Results and inputs are stored in a
`Solutions` object. Use `default_parameters` to get default model parameters. For
`MIZModel`, `init` must contain the variables `:Ei`, `:Ew`, `:h`, `:D` and `:Tg`; for
`ClassicModel`, `init` must contain `:E` and `:Tg`. A keyword argument `spectrum` must be
specified for `WIModel` to indicate the spectrum of the incident wave field.

When `lastonly=true`, only the last year of the solution is stored for each time step,
otherwise the full solution is stored. A progress bar is displayed and updated with
frequency `updatefreq`. If `updatefreq` is `Inf`, no progress bar is shown.

Refer to the documentation of the module `EnergyBalanceModel` for an example.
"""
function integrate(
    model::M, st::SpaceTime, forcing::Forcing, par::Par, init::Collection{Vec};
    lastonly::Bool=true, updatefreq::Float64=1.0
) where M<:Union{MIZModel,ClassicModel} # -> Solutions{M,F,C}
    # initialise
    vars, sols, annusol = initialise(model, st, forcing, par, init; lastonly)
    if isfinite(updatefreq)
        progress::Progress = Progress(
            length(st.T), string("Integrating ", M.name.name), updatefreq;
            infofeed=(t -> string("t = ", round(t; digits=2)))
        )
        update!(progress; feedargs=(0,))
    end # if isfinite
    # loop over time
    for ti in eachindex(st.T)
        step!(model, st.t[mod1(ti, st.nt)], forcing(st.T[ti]), vars, st, par)
        savesol!(sols, annusol, vars, ti)
        isfinite(updatefreq) && update!(progress; feedargs=(st.T[ti],))
    end # for ti
    return sols
end # function integrate

function integrate(
    model::WIModel, st::SpaceTime, forcing::Forcing, par::Par, init::Collection{Vec};
    lastonly::Bool=true, updatefreq::Float64=1.0, spectrum::Spectrum
) # -> Solutions{WIModel,F,C}
    # initialise
    if isfinite(updatefreq)
        print(SS.styled"{bold,warning:Caching wavenumber...}", lpad("", 10))
        start = time()
    end # if isfinite
    task = Threads.@spawn initialise(model, st, forcing, par, init; lastonly, spectrum)
    isfinite(updatefreq) && (
        timer = Timer(
            _ -> print("\e[10D", lpad(string(round(time()-start; digits=1)), 8), " s"),
            0.1; interval=0.1
        )
    )
    vars, sols, annusol = fetch(task)
    if isfinite(updatefreq)
        close(timer)
        println(
            "\r\e[2K",
            SS.styled"{bold,success:Wavenumber cached}",
            " in ", round(time()-start; digits=2), " s\n"
        )
        progress::Progress = Progress(
            length(st.T), "Integrating WIModel", updatefreq;
            infofeed=(t -> string("t = ", round(t; digits=2)))
        )
        update!(progress; feedargs=(0,))
    end # if isfinite
    # loop over time
    for ti in eachindex(st.T)
        step!(model, st.t[mod1(ti, st.nt)], forcing(st.T[ti]), vars, st, par; spectrum)
        savesol!(sols, annusol, vars, ti)
        isfinite(updatefreq) && update!(progress; feedargs=(st.T[ti],))
    end # for ti
    return sols
end # function integrate

"""
    solve(prob::EBMProblem; lastonly::Bool=true, updatefreq::Float64=1.0) -> Solutions{M,F,C}

Integrate the `EBMProblem` `prob` and return the results in a `Solutions` object. This is
the high-level entry point that dispatches to the appropriate `integrate` method for the
problem's model, forwarding the spectrum automatically when the model is a `WIModel`.

When `lastonly=true`, only the last year of the solution is stored for each timestep,
otherwise the full solution is stored. A progress bar is displayed and updated with
frequency `updatefreq`. If `updatefreq` is `Inf`, no progress bar is shown.

# Examples
```julia-repl
julia> problem = EBMProblem(WIModel());

julia> sols = solve(problem)
Wavenumber cached in 7.29 s

Integrating WIModel
 100000/100000 [━━━━━━━━━━━━━━━━━━━━━━━━━━━]  100%
 0:39/-0:00 2575.71/sec                     Done ✓
 t = 50.0
Solutions{WIModel, sin, false} with:
  12 solution variables: Set([:Ti, :n, :D, :h, :lambda, :phi, :Ew, :E, :Tw, :T, :Ei, :Ewave])
  on 180 latitudinal gridboxes: [0.00436331, 0.0130896 … 2, 0.999914, 0.99999]
  and 2000 timesteps: 49.00025:0.0005:49.99975
  with forcing Forcing(0.0) (constant forcing)
```
"""
solve(prob::EBMProblem; lastonly::Bool=true, updatefreq::Float64=1.0) =
    prob.model isa WIModel ?
        integrate(
            prob.model, prob.st, prob.forcing, prob.parameters, prob.initconds;
            lastonly, updatefreq, spectrum=prob.spectrum
        ) :
        integrate(
            prob.model, prob.st, prob.forcing, prob.parameters, prob.initconds;
            lastonly, updatefreq
        )

end # module Infrastructure
