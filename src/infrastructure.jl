module Infrastructure # EnergyBalanceModel.

using ..Utilities

import InteractiveUtils as IU, SparseArrays as SA, Statistics as Stats, Integrals as Intgr

export AbstractModel, ClassicModel, MIZModel, WIModel
export Collection, Forcing, Par, Solutions, SpaceTime, Vec
export default_parval
export get_diffop
export duration, hemispheric_mean, ice_area
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

function Collection(args...) # -> Coolection{V}
    dict = Dict(args...)
    return Collection{eltype(dict).types[2]}(dict)
end # function Collection

Base.getproperty(coll::Collection, key::Symbol) = getindex(getfield(coll, :dict), key) # -> V
Base.setproperty!(coll::Collection, key::Symbol, val) = setindex!(getfield(coll, :dict), val, key) # -> Dict{Symbol,V}
Base.propertynames(coll::Collection)::Set{Symbol} = Set(keys(getfield(coll, :dict)))
Base.length(coll::Collection)::Int = length(getfield(coll, :dict))
Base.hash(coll::Collection, h::UInt)::UInt = hash(getfield(coll, :dict), h)

function Base.show(io::IO, coll::Collection)::Nothing
    buffer = iobuffer(io)
    show(buffer, getfield(coll, :dict))
    str = replace(String(take!(buffer.io)), "Dict"=>string(typeof(coll)))
    print(io, str)
    return nothing
end # function Base.show

function Base.show(io::IO, ::MIME"text/plain", coll::Collection)::Nothing
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
    nt::Int # number of timesteps per year (limited by numerical stability)
    dt::Float64 # timestep
    t::Vec # time vector in a year
    winter::@NamedTuple{t::Float64, inx::Int}
    summer::@NamedTuple{t::Float64, inx::Int}

    function SpaceTime{F}(
        urange::NTuple{2,Float64}, nx::Int, nt::Int;
        winter::Float64=0.0, summer::Float64=0.5
    ) where F
        du = (urange[2]-urange[1]) / nx
        u = range(urange[1] + du/2, urange[2] - du/2, nx)
        x = F.(u)
        dt = 1 / nt
        t = collect(range(dt/2, 1 - dt/2, nt))
        winterinx = clamp(round(Int, nt*winter), 1, nt)
        summerinx = clamp(round(Int, nt*summer), 1, nt)
        return new{F}(
            nx, u, x, nt, dt, t, (t=winter, inx=winterinx), (t=summer, inx=summerinx)
        )
    end # function SpaceTime{F}
end # struct SpaceTime{F}

SpaceTime{identity}(nx::Int, nt::Int; kwargs...) = SpaceTime{identity}((0.0, 1.0), nx, nt; kwargs...)
SpaceTime{sin}(nx::Int, nt::Int; kwargs...) = SpaceTime{sin}((0.0, pi/2), nx, nt; kwargs...)
SpaceTime(args...; kwargs...) = SpaceTime{identity}(args...; kwargs...)

Base.show(io::IO, st::SpaceTime)::Nothing = print(
    io, typeof(st), '(', st.nx, ", ", st.nt, ", ", ')'
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
struct Forcing{V}
    base::Float64 # base forcing
    peak::Float64 # peak forcing
    cool::Float64 # forcing after cooldown
    rate::Float64 # rates of change
    tol::Float64 # tolerance for change
    stable_wait::Int # minimum number of years to wait for stablising before ramping up or down
    f::Vec # log of forcing in each year
    stages::Vector{Int} # ['1', 2:start_warming, 3:end_warming, 4:start_cooling, 5:end_cooling, 6:end]

    # constant forcing
    Forcing(base::Float64, tol::Float64=0.1) = new{false}(
        base, base, base, NaN, tol, -1, Vec(), [1; fill(-1, 5)]
        #     ^peak ^cool ^rate     ^stable_wait   ^stages
    )
    function Forcing(
        base::Float64, peak::Float64, cool::Float64, rate::Float64=0.2, tol::Float64=1e-3, stable_wait::Int=50
    )
        base+rate < peak > cool+rate || # check input
            throw(ArgumentError("peak is not the largest value or rate can not resolve the change."))
        new{true}(base, peak, cool, rate, tol, stable_wait, Vec(), [1; fill(-1, 5)])
    end # function Forcing
end # struct Forcing{F}

Base.getindex(f::Forcing, year::Int)::Float64 = f.f[year]

function Base.show(io::IO, forcing::Forcing{false})::Nothing
    print(io, "Forcing(", forcing.base, ')')
    printstyled(io, " (constant forcing)", color=:light_black)
    return nothing
end # function Base.show

Base.show(io::IO, forcing::Forcing{true})::Nothing = print(
    io,
    "Forcing(", forcing.base, ", ", forcing.peak, ", ", forcing.cool, ')'
)

@enum ClimateChangeState begin
    CCS_StablisingFirst = 1
    CCS_Warming = -2
    CCS_StablisingPeak = 3
    CCS_Cooling = -4
    CCS_StablisingLast = 5
    CCS_End = 6
end # enum ClimateChangeState

mutable struct ClimateChange{V}
    forcing::Forcing{V}
    year::Int
    f::Float64
    lim::Int
    avgT::Float64
    history::NTuple{2,Float64} # (previous, last)
end # struct ClimateChange

function ClimateChange!(forcing::Forcing{V}, lim::Int=(V ? 10000 : 100)) where V # TODO lim when V=true
    lim < 1 && throw(ArgumentError("Duration must be a positive integer."))
    empty!(forcing.f)
    forcing.stages = [1; fill(-1, 5)]
    return ClimateChange{V}(forcing, -1, NaN, lim, NaN, (NaN, NaN))
    #                                ^year ^f      ^avgT ^history
end # function ClimateChange!

function maintainf(cc::ClimateChange, state::ClimateChangeState)::Bool
    abs(cc.history[2]-cc.history[1]) >= cc.forcing.tol && any(>=(cc.forcing.tol)∘abs, cc.avgT .- (cc.history)) &&
        return true # unstable
    if Integer(state) > 0 # stablising states
        since = cc.forcing.stages[Integer(state)]
        cc.year - since < cc.forcing.stable_wait && return true # didn't wait for long enough
    end # if <
    return false
end # function maintainf

function next_fs!(cc::ClimateChange{true}, ::Val{CCS_StablisingFirst})::Tuple{Float64,ClimateChangeState}
    cc.forcing.stages[abs(Integer(CCS_Warming))] = cc.year + 1
    return (cc.forcing.base + cc.forcing.rate, CCS_Warming)
end # function next_fs

function next_fs!(cc::ClimateChange{true}, ::Val{CCS_Warming})::Tuple{Float64,ClimateChangeState}
    peaked = cc.forcing.peak - cc.f <= cc.forcing.rate
    if peaked
        f = cc.forcing.peak
        newstate = CCS_StablisingPeak
        cc.forcing.stages[Integer(newstate)] = cc.year
    else # !peaked
        f = cc.f + cc.forcing.rate
        newstate = CCS_Warming
    end # peaked, else
    return (f, newstate)
end # function next_fs

function next_fs!(cc::ClimateChange{true}, ::Val{CCS_StablisingPeak})::Tuple{Float64,ClimateChangeState}
    cc.forcing.stages[abs(Integer(CCS_Cooling))] = cc.year + 1
    return (cc.forcing.peak - cc.forcing.rate, CCS_Cooling)
end # function next_fs

function next_fs!(cc::ClimateChange{true}, ::Val{CCS_Cooling})::Tuple{Float64,ClimateChangeState}
    cooled = cc.f - cc.forcing.cool <= cc.forcing.rate
    if cooled
        f = cc.forcing.cool
        newstate = CCS_StablisingLast
        cc.forcing.stages[Integer(newstate)] = cc.year
    else # !cooled
        f = cc.f - cc.forcing.rate
        newstate = CCS_Cooling
    end # cooled, else
    return (f, newstate)
end # function next_fs

function Base.iterate(cc::ClimateChange)::Tuple{Float64,ClimateChangeState}
    cc.year = 1
    cc.avgT = NaN
    cc.history = (NaN, NaN)
    cc.forcing.stages[Integer(CCS_StablisingFirst)] = 1
    return (cc.forcing.base, CCS_StablisingFirst)
end

function Base.iterate(cc::ClimateChange{V}, state::ClimateChangeState) where V # -> Tuple{Float64,ClimateChangeState}/Nothing
    isnan(cc.avgT) && throw(ArgumentError("ClimateChange.avgT must be updated before iterating."))
    cc.year + 1 > cc.lim && return nothing # exceed the upper integral limit
    if maintainf(cc, state) # forcing and state do not change
        f = cc.f
        newstate = state
    elseif !V || state === CCS_StablisingLast # end of integration
        cc.forcing.stages[Integer(CCS_End)] = cc.year
        return nothing
    else # force change
        f, newstate = next_fs!(cc, Val(state))
    end # if &&, elseif ||, else
    cc.year += 1
    cc.f = f
    push!(cc.forcing.f, f)
    cc.history = (cc.history[2], cc.avgT)
    cc.avgT = NaN
    return (f, newstate)
end # function Base.iterate

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
    raw::Collection{Vector{Vec}} # solution storage
    annual::@NamedTuple{
        winter::Collection{Vector{Vec}}, summer::Collection{Vector{Vec}}, avg::Collection{Vector{Vec}}
    } # seasonal peak and annual avg

    function Solutions{M}(
        st::SpaceTime{F}, forcing::Forcing{V}, par::Par, init::Collection{Vec}, vars::Set{Symbol}
    ) where {M<:AbstractModel, F, V} # Solutions
        solraw = Collection{Vector{Vec}}() # raw solution storage
        seasonaltemp = Collection{Vector{Vec}}() # seasonal solution storage
        for var in vars
            setproperty!(solraw, var, Vector{Vec}())
            sizehint!(getproperty(solraw, var), st.nt)
            setproperty!(seasonaltemp, var, Vector{Vec}())
            sizehint!(getproperty(seasonaltemp, var), 100)
        end # for var
        return new{M,F,V}(
            st, # spacetime
            Vec(),
            forcing,
            par, # parameters
            init, # initconds
            solraw, # raw
            (
                winter=deepcopy(seasonaltemp),
                summer=deepcopy(seasonaltemp),
                avg=deepcopy(seasonaltemp)
            ) # ( # seasonal
        ) # new
    end # function Solutions
end # struct Solutions{M,F,V}

duration(sols::Solutions)::Float64 = length(sols.forcing.f)

Base.show(io::IO, sols::Solutions)::Nothing = print(
    io, typeof(sols), " with ", length(sols.raw), " variables"
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
    print(io, " for ", duration(sols), " years")
    return nothing
end # function Base.show

# default parameter values
const default_parval = Par(
    :D => 0.6, # diffusivity for heat transport (W m^-2 K^-1)
    :A => 193.0, # OLR when T = T_m (W m^-2)
    :B => 2.1, # OLR temperature dependence (W m^-2 K^-1)
    :cw => 9.8, # ocean mixed layer heat capacity (W yr m^-2 K^-1)
    :S0 => 420.0, # insolation at equator  (W m^-2)
    :S1 => 338.0, # insolation seasonal dependence (W m^-2)
    :S2 => 240.0, # insolation spatial dependence (W m^-2)
    :a0 => 0.7, # ice-free co-albedo at equator
    :a2 => 0.1, # ice-free co-albedo spatial dependence
    :ai => 0.4, # co-albedo where there is sea ice
    :Fb => 4.0, # heat flux from ocean below (W m^-2)
    :k => 2.0, # sea ice thermal conductivity (W m^-2 K^-1)
    :Lf => 9.5, # sea ice latent heat of fusion (W yr m^-3)
    :cg => 1e-3 * 9.8, # ghost layer heat capacity(W yr m^-2 K^-1)
    :tau => 1e-5 * 9.8, # ghost layer coupling timescale (yr)
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
        vecvec = getproperty(annusol.raw, var)
        length(vecvec) == annusol.spacetime.nt ||
            throw(
                ArgumentError(
                    "Length of raw solution vector for $var does not match the number of timesteps per year, when calculating annual mean."
                )
            ) # throw
        setproperty!(means, var, crossmean(vecvec))
    end # for var
    return means
end # function annual_mean

function save_step!(annusol::Solutions, vars::Collection{Vec}) # -> Solutions{M,F,V}
    varscp = deepcopy(vars) # avoid reference issues
    # save raw data to annual
    foreach(
        var -> push!(getproperty(annusol.raw, var), getproperty(varscp, var)), # !
        propertynames(annusol.raw)
    )
    return annusol
end # function save_step!

function save_year!(
    sols::Solutions, annusol::Solutions, year::Int, saveraw::Bool=false
)
    if saveraw # save raw data
        foreach(
            var -> append!(getproperty(sols.raw, var), getproperty(annusol.raw, var)),
            propertynames(sols.raw)
        )
        append!(sols.ts, year .+ sols.spacetime.t)
    end # if saveraw
    # save seasonal data
    for season in (:winter, :summer)
        foreach(propertynames(getproperty(sols.annual, season))) do var
            push!(
                getproperty(getproperty(sols.annual, season), var),
                getproperty(annusol.raw, var)[getproperty(sols.spacetime, season).inx]
            )
        end # do var
    end # for season
    # save annual mean data
    means = annual_mean(annusol)
    foreach(
        var -> push!(getproperty(sols.annual.avg, var), getproperty(means, var)),
        propertynames(sols.annual.avg)
    )
    # empty annual raw data
    foreach(var -> empty!(getproperty(annusol.raw, var)), propertynames(annusol.raw))
    return sols
end # function save_year!


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

# TODO doc string
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

ice_area(sols::Solutions{ClassicModel}, season::Symbol, year::Int)::Float64 =
    ice_area((getproperty(sols.annual, season).E[year].<0), sols.spacetime.x)
ice_area(sols::Solutions{<:Union{MIZModel,WIModel}}, season::Symbol, year::Int)::Float64 =
    ice_area(getproperty(sols.annual, season).phi[year], sols.spacetime.x)

# stub for functions for each model
function step! end
function initialise end

"""
    integrate(model<:AbstractModel, st::SpaceTime, forcing::Forcing{false}, par::Par, init::Collection{Vec}; lastonly::Bool=true, progress::Bool=true, spectrum::Union{Spectrum,Nothing}=nothing) -> Solutions{M,_,false}

Integrate the specified model over the given `SpaceTime` with climate `Forcing`, model
parameters `par`, and initial conditions `init`. Results and inputs are stored in a
`Solutions` object. Use `default_parameters` to get default model parameters. For
`MIZModel`, `init` must contain the variables `:Ei`, `:Ew`, `:h`, `:D` and `:Tg`; for
`ClassicModel`, `init` must contain `:E` and `:Tg`. A keyword argument `spectrum` must be
specified for `WIModel` to indicate the spectrum of the incident wave field.

When `lastonly=true`, only the last year of the solution is stored for each time step,
otherwise the full solution is stored. A progress bar is displayed if `progress=true`.

Refer to the documentation of the module `EnergyBalanceModel` for an example.
""" # TODO update doc string
function integrate!(
    model::M, st::SpaceTime, forcing::Forcing, par::Par, init::Collection{Vec}, lim::Union{<:Integer,Nothing}=nothing;
    lastonly::Bool=true, progress::Bool=true, spectrum::S=nothing
) where {M<:AbstractModel, S} # -> Solutions{M,F,false}
    # check the type of spectrum
    S === Nothing || Symbol(S) === :Spectrum ||
        throw(ArgumentError("Keyword argument `spectrum` must be of type Spectrum or Nothing."))
    # check lim
    isnothing(lim) || lim > 0 ||
        throw(ArgumentError("lim must be a positive integer or nothing."))
    # warn if spectrum is provided for non-WIModel
    (M === WIModel || S===Nothing) ||
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
end # function integrate!

function integrate(
    model::AbstractModel, st::SpaceTime, forcing::Forcing, par::Par, init::Collection{Vec},
    lim::Union{<:Integer,Nothing}=nothing;
    lastonly::Bool=true, progress::Bool=true, spectrum=nothing
) # -> Solutions{M,F,true}
    dcpf = deepcopy(forcing)
    return integrate!(model, st, dcpf, par, init, lim; lastonly, progress, spectrum)
end # function integrate

end # module Infrastructure
