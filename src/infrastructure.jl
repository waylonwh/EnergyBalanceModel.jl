module Infrastructure # EnergyBalanceModel.

using ..Utilities

import EnergyBalanceModel
import SparseArrays as SA, Statistics as Stats

export Vec, Collection, SpaceTime, Solutions, Forcing
export default_parval, miz_paramset, classic_paramset
export default_parameters, get_diffop, diffusion!, D∇²!, diffusion, D∇², annual_mean
export integrate

const Vec = Vector{Float64} # abbreviation for vector type used in model

struct Collection{V}
    dict::Dict{Symbol,V}
    Collection{V}(args...) where V = new(Dict{Symbol,V}(args...))
end # struct Collection

(Base.getproperty(coll::Collection{V}, key::Symbol)::V) where V = getindex(getfield(coll, :dict), key)
(Base.setproperty!(coll::Collection{V}, key::Symbol, val::V)::Dict{Symbol,V}) where V =
    setindex!(getfield(coll, :dict), val, key)
(Base.keys(coll::Collection{V})::Base.KeySet{Symbol, Dict{Symbol,V}}) where V = keys(getfield(coll, :dict))
(Base.length(coll::Collection{V})::Int) where V = length(getfield(coll, :dict))

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

struct SpaceTime{F}
    nx::Int # number of evenly spaced latitudinal gridboxes (equator to pole)
    u::Vec # grid before modification/scale
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
        dx = (urange[2]-urange[1]) / nx
        u = collect(urange[1] + dx/2.0 : dx : urange[2] - dx/2.0)
        x = F.(u)
        dt = 1.0 / nt
        t = collect(range(dt/2.0, 1.0 - dt/2.0, nt))
        T = dt/2.0 : dt : dur - dt/2.0
        winterinx = round(Int, nt*winter)
        summerinx = round(Int, nt*summer)
        return new{F}(
            nx, u, x, dur, nt, dt, t, T, (t=winter, inx=winterinx), (t=summer, inx=summerinx)
        )
    end # function SpaceTime{F}
end # struct SpaceTime{F}

SpaceTime{identity}(nx::Int, nt::Int, dur::Int; kwargs...) = SpaceTime{identity}((0.0, 1.0), nx, nt, dur; kwargs...)
SpaceTime{sin}(nx::Int, nt::Int, dur::Int; kwargs...) = SpaceTime{sin}((0.0, pi/2.0), nx, nt, dur; kwargs...)
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

struct Forcing{F}
    base::Float64 # base forcing
    peak::Float64 # peak forcing
    cool::Float64 # forcing after cooldown
    holdyrs::NTuple{2,Int} # years to hold at (base, peak) forcing
    rates::NTuple{2,Float64} # rates of change
    domain::NTuple{5,Int} # years at which forcing pattern changes

    # constant forcing
    Forcing(base::Float64) = new{true}(
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
        warming = (peak - base) / rates[1]
        rates[1]>0 && isinteger(warming) ?
            @.(domainvec[3:5] += warming) :
            throw(ArgumentError("Warming time must be positive integer. Got $warming y."))
        # hold at peak
        @. domainvec[4:5] += holdyrs[2]
        # time to cool
        cooling = (cool - peak) / rates[2]
        rates[2]<0 && isinteger(cooling) ?
            domainvec[5] += cooling :
            throw(ArgumentError("Cooling time must be positive integer. Got $cooling y."))
        return new{false}(base, peak, cool, holdyrs, rates, Tuple(domainvec))
    end # function Forcing
end # struct Forcing{F}

function Base.show(io::IO, forcing::Forcing{true})::Nothing
    print(io, typeof(forcing), '(', forcing.base, ')')
    printstyled(io, " (constant forcing)", color=:light_black)
    return nothing
end # function Base.show

Base.show(io::IO, forcing::Forcing{false})::Nothing = print(
    io,
    typeof(forcing), '(', forcing.base, " ↗ ", forcing.peak, " ↘ ", forcing.cool, ')'
)

function Base.show(io::IO, ::MIME"text/plain", forcing::Forcing{true})::Nothing
    println(io, typeof(forcing), '(', forcing.base, ") is constant:")
    print(io, "  F(t)=", forcing.base, ", t∈[0,∞)")
    return nothing
end # function Base.show

function Base.show(io::IO, ::MIME"text/plain", forcing::Forcing{false})::Nothing
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
        ' ', rate>0.0 ? '+' : '-', ' ',
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
(forcing::Forcing{true})(::Float64)::Float64 = forcing.base # constant forcing
function (forcing::Forcing{false})(T::Float64)::Float64 # varying forcing
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
end # function (forcing::Forcing{false})

struct Solutions{F,C}
    spacetime::SpaceTime{F} # space and time which solutions are defined on
    ts::Vec # time vector for stored solution
    forcing::Forcing{C} # climate forcing
    parameters::Collection{Float64} # model parameters
    initconds::Collection{Vec} # initial conditions
    lastonly::Bool # store only last year of solution
    debug::Union{Expr,Nothing} # store debug variables
    raw::Collection{Vector{Vec}} # solution storage
    seasonal::@NamedTuple{
        winter::Collection{Vector{Vec}}, summer::Collection{Vector{Vec}}, avg::Collection{Vector{Vec}}
    } # seasonal peak and annual avg

    function Solutions(
        st::SpaceTime{F}, forcing::Forcing{C}, par::Collection{Float64}, init::Collection{Vec}, vars::Set{Symbol},
        lastonly::Bool=true;
        debug::Union{Expr,Nothing}=nothing
    ) where {F, C} # Solutions
        if lastonly
            dur_store = 1
            ts::Vec = st.dur-1.0 + st.dt/2.0 : st.dt : st.dur - st.dt/2.0
        else # !lastonly
            dur_store = st.dur
            ts = st.dt/2.0 : st.dt : st.dur - st.dt/2.0
        end # if lastonly, else
        if !isnothing(debug)
            push!(vars, :debug)
        end # if !=
        # construct raw solution storage
        solraw = Collection{Vector{Vec}}()
        foreach((var -> setproperty!(solraw, var, Vector{Vec}(undef, length(ts)))), vars)
        # construct seasonal solution storage template
        seasonaltemp = Collection{Vector{Vec}}()
        foreach((var -> setproperty!(seasonaltemp, var, Vector{Vec}(undef, st.dur))), vars)
        return new{F, C}(
            st, # spacetime
            ts,
            forcing,
            par, # parameters
            init, # initconds
            lastonly,
            debug,
            solraw, # raw
            (
                winter=deepcopy(seasonaltemp),
                summer=deepcopy(seasonaltemp),
                avg=deepcopy(seasonaltemp)
            ) # ( # seasonal
        ) # new
    end # function Solutions
end # struct Solutions{F,C}

(Base.show(io::IO, sols::Solutions{F,C})::Nothing) where {F,C} = print(
    io,
    typeof(sols), '(',
    sols.spacetime.nx, '×', length(sols.ts), "@(", first(sols.ts), ':', sols.spacetime.dt, ':', last(sols.ts), "), ",
    keys(sols.raw),
    ')'
)

function Base.show(io::IO, ::MIME"text/plain", sols::Solutions{F,C})::Nothing where {F,C}
    println(io, typeof(sols), " with:")
    println(io, "  ", length(sols.raw), " solution variables: ", keys(sols.raw))
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
const default_parval = Collection{Float64}(
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
    :F => 0.0, # radiative forcing (W m^-2)
    :cg => 0.01 * 9.8, # ghost layer heat capacity(W yr m^-2 K^-1)
    :tau => 1e-5, # ghost layer coupling timescale (yr)
    :Tm => 0.0, # mean temperature (C)
    :m1 => 1.6e-6 * 31536000, # empirical constants of lateral melt
    :m2 => 1.36, # empirical constants of lateral melt
    :alpha => 0.66, # floe geometry constant, Ai = alpha * D^2
    :rl => 0.5, # lead region width (m)
    :Dmin => 1.0, # new pancake size (m)
    :Dmax => 156, # largest floe length (m)
    :hmin => 0.1, # new pancake thickness (m)
    :kappa => 0.01 * 31536000 # floe welding parameter
) # Collection{Float64}

# parameters used in each model
const miz_paramset = Set{Symbol}(
    (
        :D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :Tm, :m1, :m2, :alpha,
        :rl, :Dmin, :Dmax, :hmin, :kappa
    )
)
const classic_paramset = Set{Symbol}(
    (:D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :F, :cg, :tau)
)

# Create a parameter dictionary from default values for a given Set
function default_parameters(paramset::Set{Symbol})::Collection{Float64}
    setvec = collect(paramset)
    return Collection{Float64}(setvec .=> getproperty.(Ref(default_parval), setvec))
end # function get_defaultpar
default_parameters(model::Symbol)::Collection{Float64} =
    model===:MIZ ? default_parameters(miz_paramset) : default_parameters(classic_paramset)

# calculate diffusion operator matrix
@persistent(
    diffop::SA.SparseMatrixCSC{Float64,Int64} = SA.spzeros(Float64, 0, 0),

    @inline function get_diffop(nx::Int)::SA.SparseMatrixCSC{Float64,Int64}
        if size(diffop) != (nx, nx) # recalculate diffusion operator
            dx = 1.0 / nx
            xb = dx:dx:1.0-dx
            lambda = @. (1 - xb^2) / dx^2
            l1 = pushfirst!(-copy(lambda), 0.0)
            l2 = push!(-copy(lambda), 0.0)
            l3 = -l1 - l2
            diffop = SA.spdiagm(-1 => -l1[2:nx], 0 => -l3, 1 => -l2[1:nx-1])
        end # if !=
        return diffop
    end # function get_diffop
) # @persistent

# diffusion for equal spaced grid
@inline (diffusion!(
    base::Vector{T}, temp::Vector{T}, st::SpaceTime{identity}, par::Collection{Float64}
)::Vector{T}) where T<:Number = base .+= par.D * get_diffop(st.nx) * temp

# diffusion for non-equal spaced grid
@persistent(
    diffx::Vector{Float64}, mxxph::Vector{Float64}, mxxmh::Vector{Float64},
    phmmh::Vector{Float64}, i::UnitRange{Int},
    xid::UInt64=UInt64(0),

    @inline function diffusion!(
        base::Vector{T}, temp::Vector{T}, st::SpaceTime{F}, par::Collection{Float64}
    )::Vector{T} where {T<:Number,F}
        # store x if changed
        if xid != objectid(st.x)
            x = [-st.x[1]; st.x; 2-st.x[end]]
            diffx = diff(x)
            diffT = zeros(Float64, st.nx+1)
            i = 2:st.nx+1
            xxph = @. (x[i+1]+x[i]) / 2.0
            xxmh = @. (x[i]+x[i-1]) / 2.0
            mxxph = @. 1.0 - xxph^2
            mxxmh = @. 1.0 - xxmh^2
            phmmh = @. xxph - xxmh
            xid = objectid(st.x)
        end # if !=
        diffT = Vector{T}(undef, st.nx+1) # TODO eliminate memory allocations?
        @inbounds diffT[1] = diffT[end] = zero(T)
        @inbounds diffT[2:st.nx] .= diff(temp) # !
        @inbounds @. base += par.D * (mxxph * diffT[i]/diffx[i] - mxxmh * diffT[i-1]/diffx[i-1]) / phmmh # !
        return base
    end # function diffusion!
) # @persistent

@inline (diffusion(T::Vec, st::SpaceTime{F}, par::Collection{Float64})::Vec) where F =
    diffusion!(zeros(Float64, length(T)), T, st, par)

const D∇² = diffusion
const D∇²! = diffusion!

# calculate annual mean
function annual_mean(annusol::Solutions{F,C})::Collection{Vec} where {F, C}
    # calculate annual mean for each variable except temperatures
    means = Collection{Vec}()
    foreach(
        (var -> setproperty!(means, var, crossmean(getproperty(annusol.raw, var)))),
        keys(annusol.raw)
    )
    return means
end # function annual_mean

(annual_mean(forcing::Forcing{C}, st::SpaceTime{F}, year::Int)::Float64) where {C, F} =
    Stats.mean(forcing.(year-1 .+ st.t))

function savesol!(
    sols::Solutions{F,C}, annusol::Solutions{F,C}, vars::Collection{Vec}, tinx::Int
)::Solutions{F,C} where {F, C}
    varscp = deepcopy(vars) # avoid reference issues
    year = ceil(Int, sols.spacetime.T[tinx])
    ti = mod1(tinx, sols.spacetime.nt) # index of time in the year
    # save raw data to annual
    foreach(
        (var -> getproperty(annusol.raw, var)[ti] = getproperty(varscp, var)), # !
        keys(annusol.raw)
    )
    # save raw data
    if !sols.lastonly # save all raw data
        foreach(
            (var -> setindex!(getproperty(sols.raw, var), getproperty(varscp, var), tinx)),
            keys(sols.raw)
        )
    elseif tinx > length(sols.spacetime.T) - sols.spacetime.nt # save the raw data of the last year
        foreach(
            (var -> setindex!(getproperty(sols.raw, var), getproperty(varscp, var), ti)),
            keys(sols.raw)
        )
    end # if !, elseif
    # save seasonal data
    if ti == sols.spacetime.winter.inx
        foreach(
            (var -> setindex!(getproperty(sols.seasonal.winter, var), getproperty(varscp, var), year)),
            keys(sols.seasonal.winter)
        )
    elseif ti == sols.spacetime.summer.inx
        foreach(
            (var -> setindex!(getproperty(sols.seasonal.summer, var), getproperty(varscp, var), year)),
            keys(sols.seasonal.summer)
        )
    elseif ti == sols.spacetime.nt # calculate annual average
        means = annual_mean(annusol)
        foreach(
            (var -> setindex!(getproperty(sols.seasonal.avg, var), getproperty(means, var), year)),
            keys(sols.seasonal.avg)
        )
    end # if ==, elseif*2
    return sols
end # function savesol!

function integrate(
    model::Symbol, st::SpaceTime{F}, forcing::Forcing{C}, par::Collection{Float64}, init::Collection{Vec};
    lastonly::Bool=true, debug::Union{Expr,Nothing}=nothing, verbose::Bool=false
)::Solutions{F,C} where {F, C}
    # initialise
    vars = deepcopy(init)
    solvars = Set{Symbol}((:E, :T, :h)) # always solve for these
    if model === :MIZ # add MIZ variables
        union!(solvars, Set{Symbol}((:Ei, :Ew, :Ti, :Tw, :D, :phi, :n)))
    end # if ===
    Modu::Module = EnergyBalanceModel.eval(model)
    sols = Solutions(st, forcing, par, init, solvars, lastonly; debug=debug)
    annusol = Solutions(st, forcing, par, init, solvars, true; debug=debug) # for calculating annual means
    progress = Progress(length(st.T), "Integrating"; infofeed=(t -> string("t = ", round(t, digits=2))))
    update!(progress; feedargs=(0,))
    # loop over time
    for ti in eachindex(st.T)
        Modu.step!(st.t[mod1(ti, st.nt)], forcing(st.T[ti]), vars, st, par; debug=debug, verbose=verbose)
        savesol!(sols, annusol, vars, ti)
        update!(progress; feedargs=(st.T[ti],))
    end # for ti
    return sols
end # function integrate

end # module Infrastructure
