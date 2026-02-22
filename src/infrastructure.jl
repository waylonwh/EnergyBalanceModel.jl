module Infrastructure # EnergyBalanceModel.

using ..Utilities

import InteractiveUtils as IU, SparseArrays as SA, Statistics as Stats

export AbstractModel, ClassicModel, MIZModel, WIModel
export Collection, Forcing, Par, Solutions, SpaceTime, Vec
export classic_paramset, default_parameters, default_parval, miz_paramset
export get_diffop
export annual_mean, hemispheric_mean
export integrate

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

julia> parameters.F = 0.0; parameters.F
0.0
```
"""
struct Collection{V}
    dict::Dict{Symbol,V}
    Collection{V}(args...) where V = new(Dict{Symbol,V}(args...))
end # struct Collection

(Base.getproperty(coll::Collection{V}, key::Symbol)::V) where V = getindex(getfield(coll, :dict), key)
(Base.setproperty!(coll::Collection{V}, key::Symbol, val::V)::Dict{Symbol,V}) where V =
    setindex!(getfield(coll, :dict), val, key)
(Base.propertynames(coll::Collection{V})::Set{Symbol}) where V = Set(keys(getfield(coll, :dict)))
(Base.length(coll::Collection{V})::Int) where V = length(getfield(coll, :dict))
(Base.hash(coll::Collection{V}, h::UInt)::UInt) where V = hash(getfield(coll, :dict), h)

function Base.show(io::IO, coll::Collection{V})::Nothing where V
    buffer = iobuffer(io)
    show(buffer, getfield(coll, :dict))
    str = replace(String(take!(buffer.io)), "Dict"=>string(typeof(coll)))
    print(io, str)
    return nothing
end # function Base.show

function Base.show(io::IO, ::MIME"text/plain", coll::Collection{V})::Nothing where V
    buffer = iobuffer(io)
    show(buffer, MIME("text/plain"), getfield(coll, :dict))
    str = replace(
        String(take!(buffer.io)),
        string(typeof(getfield(coll, :dict))) => string(typeof(coll))
    )
    print(io, str)
    return nothing
end # function Base.show

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
    dur::Int # duration of simulation in years # TODO
    nt::Int # number of timesteps per year (limited by numerical stability)
    dt::Float64 # timestep
    t::Vec # time vector in a year
    T::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64} # full time series
    winter::@NamedTuple{t::Float64, inx::Int}
    summer::@NamedTuple{t::Float64, inx::Int}

    function SpaceTime{F}(
        urange::NTuple{2,Float64}, nx::Int, nt::Int, dur::Int;
        winter::Float64=0.26125, summer::Float64=0.77375 # TODO definition of "winter" and "summer"
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

(Base.show(io::IO, st::SpaceTime{F})::Nothing) where F = print(
    io,
    typeof(st), '(', st.nx, ", ", st.nt, ", ", st.dur, ')'
)

function Base.show(io::IO, ::MIME"text/plain", st::SpaceTime{F})::Nothing where F
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

    Forcing(base::Float64, peak::Float64, cool::Float64, step::Float64=0.1, tol::Float64)

Defines a time varying self adaptive climate forcing that first reaches a steady state
under the `base` forcing. This takes at least 10 model years. It then gradually ramps up to
`peak` at a rate no greater than `step` per year. A step is applied only if the change in
the annual hemispheric mean temperature is less than `tol`, so the system is assured to be
close to steady state before warming. After reaching `peak`, the forcing is maintained at
`peak` for 10 years. The forcing then ramps down to `cool`, following the same step and
tolerance criteria.

`Forcing.f` stores the forcing value for each year of the simulation, and `Forcing.tip` is
the year when the forcing pattern changes from ramping up to ramping down.

# Examples
```julia-repl

```
""" # TODO Examples for Forcing
mutable struct Forcing{C}
    base::Float64 # base forcing
    peak::Float64 # peak forcing
    cool::Float64 # forcing after cooldown
    rate::Float64 # rates of change
    tol::Float64 # tolerance for change
    f::Vec # log of forcing in each year
    tip::Int # the year when the forcing pattern changes

    # constant forcing
    Forcing(base::Float64) = new{true}(
        base, base, base, NaN, NaN, Vec(), 0, 0
    )
    # warming/cooling forcing # TODO determine tol
    Forcing(base::Float64, peak::Float64, cool::Float64, rate::Float64=0.1, tol::Float64=0.1) =
        new{false}(base, peak, cool, rate, tol, Vec(), 0, 0)
end # struct Forcing{F}

function Base.show(io::IO, forcing::Forcing{true})::Nothing
    print(io, "Forcing(", forcing.base, ')')
    printstyled(io, " (constant forcing)", color=:light_black)
    return nothing
end # function Base.show

Base.show(io::IO, forcing::Forcing{false})::Nothing = print(
    io,
    "Forcing(", forcing.base, ", ", forcing.peak, ", ", forcing.cool, ')'
)

@enum ClimateChangeStage CCS_Stablising CCS_Warming CCS_Waiting CCS_Cooling CCS_End

struct ClimateChangeState{S}
    year::Int
    stable::Bool
end

mutable struct ClimateChange
    forcing::Forcing{false}
    f::Float64
    avgT::Float64
    history::NTuple{2,Float64} # (previous, last)

    ClimateChange(forcing::Forcing{false}) = new(forcing, forcing.base, NaN, (NaN, NaN))
end # struct ClimateChange

Base.iterate(cc::ClimateChange)::Tuple{Float64,ClimateChangeState{CCS_Stablising}} = (
    cc.f, ClimateChangeState{CCS_Stablising}(0, false)
)

function Base.iterate(
    cc::ClimateChange, state::ClimateChangeState{S}
)::Tuple{Float64,ClimateChangeState{T}} where {S, T}
    isnan(cc.avgT) && throw(ArgumentError("ClimateChange.avgT must be updated before iterating."))
    # TODO pick up from here
end # function Base.iterate

Base.isdone(::ClimateChange, state::ClimateChangeState) = state.stage === CCS_End

"""
    Solutions{M,F,C}

An object to store model solutions. Type parameter `M` is the model type (`MIZModel` or
`ClassicModel`); `F` is the function used to map the uniform grid to the model grid in
`SpaceTime{F}`; `C` is a boolean indicating whether the climate forcing is constant.
`C` is `true` for constant forcing.

# Fields
- `model::M`: model type
- `spacetime::SpaceTime{F}`: space and time on which solutions are defined
- `ts::Vec`: time vector for stored solutions
- `forcing::Forcing{C}`: climate forcing
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
struct Solutions{M<:AbstractModel,F,C}
    spacetime::SpaceTime{F} # space and time which solutions are defined on
    ts::Vec # time vector for stored solution
    forcing::Forcing{C} # climate forcing
    parameters::Par # model parameters
    initconds::Collection{Vec} # initial conditions
    lastonly::Bool # store only last year of solution
    raw::Collection{Vector{Vec}} # solution storage
    annual::@NamedTuple{
        winter::Collection{Vector{Vec}}, summer::Collection{Vector{Vec}}, avg::Collection{Vector{Vec}}
    } # seasonal peak and annual avg

    function Solutions{M}(
        st::SpaceTime{F}, forcing::Forcing{C}, par::Par, init::Collection{Vec},
        vars::Set{Symbol}, lastonly::Bool=true
    ) where {M<:AbstractModel, F, C} # Solutions
        if lastonly
            dur_store = 1
            ts::Vec = st.dur-1 + st.dt/2 : st.dt : st.dur - st.dt/2
        else # !lastonly
            dur_store = st.dur
            ts = st.dt/2 : st.dt : st.dur - st.dt/2
        end # if lastonly, else
        # construct raw solution storage
        solraw = Collection{Vector{Vec}}()
        foreach((var -> setproperty!(solraw, var, Vector{Vec}(undef, length(ts)))), vars)
        # construct seasonal solution storage template
        seasonaltemp = Collection{Vector{Vec}}()
        foreach((var -> setproperty!(seasonaltemp, var, Vector{Vec}(undef, st.dur))), vars)
        return new{M,F,C}(
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
            ) # ( # seasonal
        ) # new
    end # function Solutions
end # struct Solutions{M,F,C}

(Base.show(io::IO, sols::Solutions{<:AbstractModel,F,C})::Nothing) where {F, C} = print(
    io,
    typeof(sols), '(',
    sols.spacetime.nx, '×', length(sols.ts), "@(", first(sols.ts), ':', sols.spacetime.dt, ':', last(sols.ts), "), ",
    propertynames(sols.raw),
    ')'
)

function Base.show(io::IO, ::MIME"text/plain", sols::Solutions{<:AbstractModel,F,C})::Nothing where {F, C}
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

# default parameter values
let cw::Float64 = 9.8
    global const default_parval = Par(
        :D => 0.6, # diffusivity for heat transport (W m^-2 K^-1)
        :A => 193.0, # OLR when T = T_m (W m^-2)
        :B => 2.1, # OLR temperature dependence (W m^-2 K^-1)
        :cw => cw, # ocean mixed layer heat capacity (W yr m^-2 K^-1)
        :S0 => 420.0, # insolation at equator  (W m^-2)
        :S1 => 338.0, # insolation seasonal dependence (W m^-2)
        :S2 => 240.0, # insolation spatial dependence (W m^-2)
        :a0 => 0.7, # ice-free co-albedo at equator
        :a2 => 0.1, # ice-free co-albedo spatial dependence
        :ai => 0.4, # co-albedo where there is sea ice
        :Fb => 4.0, # heat flux from ocean below (W m^-2)
        :k => 2.0, # sea ice thermal conductivity (W m^-2 K^-1)
        :Lf => 9.5, # sea ice latent heat of fusion (W yr m^-3)
        :cg => 1e-3 * cw, # ghost layer heat capacity(W yr m^-2 K^-1)
        :tau => 1e-5 * cw, # ghost layer coupling timescale (yr)
        :Tm => 0.0, # mean temperature (C)
        :m1 => 1.6e-6 * 31536000, # empirical constants of lateral melt (m y^-1 K^-1)
        :m2 => 1.36, # empirical constants of lateral melt
        :alpha => 0.66, # floe geometry constant, Ai = alpha * D^2
        :rl => 0.5, # lead region width (m)
        :Dmin => 1.0, # new pancake size (m)
        :Dmax => 500.0, # largest floe length (m)
        :hmin => 0.1, # new pancake thickness (m)
        :kappa => 0.01 * 31536000, # floe welding parameter (m^2 s^-1)
        :Y => 5.5, # Effective Young's modulus (GPa)
        :nu => 0.3, # Poisson's ratio
        :rhow => 1025.0, # Water density (kg/m^3)
        :g => 9.81, # Gravitational acceleration (m/s^2),
        :Ec => 7.05e-5, # Breaking significant strain
        :Gamma => 13.0, # Viscous damping parameter (Pa m s^-1)
        :gamma => 2 + log2(0.9), # Power law exponent for floe size distribution
        :dmn => 20.0, # Chosen minmum floe diameter for the truncated power-law FSD in WIM (m)
    ) # Par
end # let cw

# parameters used in each model
const classicmodel_parvars = Set{Symbol}(
    (:D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :cg, :tau)
)
const mizmodel_parvars = push!(
    deepcopy(classicmodel_parvars),
    :Tm, :m1, :m2, :alpha, :rl, :Dmin, :Dmax, :hmin, :kappa
)
const wimodel_parvars = push!(
    deepcopy(mizmodel_parvars),
    :Y, :nu, :rhow, :g, :Ec, :Gamma, :gamma, :dmn
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
    @inline function get_diffop(st::SpaceTime{F})::SA.SparseMatrixCSC{Float64,Int64} where F
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

function annual_mean(annusol::Solutions{<:AbstractModel,F,C})::Collection{Vec} where {F, C}
    # calculate annual mean for each variable except temperatures
    means = Collection{Vec}()
    foreach(
        (var -> setproperty!(means, var, crossmean(getproperty(annusol.raw, var)))),
        propertynames(annusol.raw)
    )
    return means
end # function annual_mean

"""
    annual_mean(forcing::Forcing{C}, st::SpaceTime{F}, year::Int) -> Float64

Calculate the annual mean of the climate forcing for a given year.

# Examples
```julia-repl
julia> forcing = Forcing(0.0, 10.0, -5.0, (20, 10), (0.5, -0.5));

julia> st = SpaceTime{sin}(180, 2000, 30);

julia> annual_mean(forcing, st, 24)
1.75
```
"""
(annual_mean(forcing::Forcing{C}, st::SpaceTime{F}, year::Int)::Float64) where {C, F} =
    Stats.mean(forcing.(year-1 .+ st.t))

function savesol!(
    sols::Solutions{M,F,C}, annusol::Solutions{M,F,C}, vars::Collection{Vec}, tinx::Int
)::Solutions{M,F,C} where {M<:AbstractModel, F, C}
    varscp = deepcopy(vars) # avoid reference issues
    year = ceil(Int, sols.spacetime.T[tinx])
    ti = mod1(tinx, sols.spacetime.nt) # index of time in the year
    # save raw data to annual
    foreach(
        (var -> getproperty(annusol.raw, var)[ti] = getproperty(varscp, var)), # !
        propertynames(annusol.raw)
    )
    # save raw data
    if !sols.lastonly # save all raw data
        foreach(
            (var -> setindex!(getproperty(sols.raw, var), getproperty(varscp, var), tinx)),
            propertynames(sols.raw)
        )
    elseif tinx > length(sols.spacetime.T) - sols.spacetime.nt # save the raw data of the last year
        foreach(
            (var -> setindex!(getproperty(sols.raw, var), getproperty(varscp, var), ti)),
            propertynames(sols.raw)
        )
    end # if !, elseif
    # save seasonal data
    if ti == sols.spacetime.winter.inx
        foreach(
            (var -> setindex!(getproperty(sols.annual.winter, var), getproperty(varscp, var), year)),
            propertynames(sols.annual.winter)
        )
    elseif ti == sols.spacetime.summer.inx
        foreach(
            (var -> setindex!(getproperty(sols.annual.summer, var), getproperty(varscp, var), year)),
            propertynames(sols.annual.summer)
        )
    elseif ti == sols.spacetime.nt # calculate annual average
        means = annual_mean(annusol)
        foreach(
            (var -> setindex!(getproperty(sols.annual.avg, var), getproperty(means, var), year)),
            propertynames(sols.annual.avg)
        )
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
    int = zero(Float64)
    for i in 1:length(x)-1
        @inbounds int += (vec[i] + vec[i+1]) * (x[i+1] - x[i]) / 2
    end # for i
    return int
end # function hemispheric_mean

# stub for functions for each model
function step! end
function initialise end

"""
    integrate(model::M, st::SpaceTime{F}, forcing::Forcing{false}, par::Par, init::Collection{Vec}; lastonly::Bool=true, progress::Bool=true, spectrum::Union{Spectrum,Nothing}=nothing)::Solutions{M,F,C}

Integrate the specified model over the given `SpaceTime` with climate `Forcing`, model
parameters `par`, and initial conditions `init`. Results and inputs are stored in a
`Solutions` object. Use `default_parameters` to get default model parameters. For
`MIZModel`, `init` must contain the variables `:Ei`, `:Ew`, `:h`, `:D` and `:Tg`; for
`ClassicModel`, `init` must contain `:E` and `:Tg`. A keyword argument `spectrum` must be
specified for `WIModel` to indicate the spectrum of the incident wave field.

When `lastonly=true`, only the last year of the solution is stored for each time step,
otherwise the full solution is stored. A progress bar is displayed if `progress=true`.

Refer to the documentation of the module `EnergyBalanceModel` for an example.
"""
function integrate(
    model::M, st::SpaceTime{F}, forcing::Forcing{true}, par::Par, init::Collection{Vec};
    lastonly::Bool=true, progress::Bool=true,
    spectrum#=::Union{Spectrum,Nothing}=#=nothing
)::Solutions{M,F} where {M<:AbstractModel, F}
    # warn if spectrum is provided for non-WIModel
    (M === WIModel || isnothing(spectrum)) ||
        @warn "Keyword argument `spectrum` is ignored as $M does not have a WIM component."
    # initialise
    forcing.f = fill(forcing.base, st.dur) # fill in the forcing value for each year
    vars, sols, annusol = initialise(model, st, forcing, par, init; lastonly)
    prog::Progress = Progress(
        st.T[last], "t", "Integrating";
        infofeed=((t, tps) -> string("t = ", round(t; digits=2), "  ", round(1/tps; digits=2), "itr/s"))
    )
    progress && start!(prog; feedargs=(0,))
    # loop over time
    for ti in eachindex(st.T)
        vt = @timed step!(model, st.t[mod1(ti, st.nt)], forcing(st.T[ti]), vars, st, par; spectrum)
        savesol!(sols, annusol, vars, ti)
        t = st.T[ti]
        progress && update!(prog, t; feedargs=(t, vt.time))
    end # for ti
    return sols
end # function integrate

function integrate(
    model::M, st::SpaceTime{F}, forcing::Forcing{false}, par::Par, init::Collection{Vec};
    lastonly::Bool=true, progress::Bool=true,
    spectrum#=::Union{Spectrum,Nothing}=#=nothing
)::Solutions{M,F} where {M<:AbstractModel, F} # constant forcing
    # warn if spectrum is provided for non-WIModel
    (M === WIModel || isnothing(spectrum)) ||
        @warn "Keyword argument `spectrum` is ignored as $M does not have a WIM component."
    # initialise
    vars, sols, annusol = initialise(model, st, forcing, par, init; lastonly)
    # reach base equilibrium
    prog::Progress = Progress(
        forcing.tol, forcing.tol, "|ΔT|", "Reaching Steady State";
        infofeed=((t, tps) -> string("t = ", round(t; digits=2), "  ", round(1/tps; digits=2), "itr/s"))
    )
    progress && start!(prog; feedargs=(0,))
    # loop over time
    for ti in eachindex(st.T)
        vt = @timed step!(model, st.t[mod1(ti, st.nt)], forcing(st.T[ti]), vars, st, par; spectrum)
        savesol!(sols, annusol, vars, ti)
        t = st.T[ti]
        progress && update!(prog, t; feedargs=(t, vt.time))
    end # for ti
    return sols
end # function integrate

end # module Infrastructure
