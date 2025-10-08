module EnergyBalanceModel # .


module Utilities # EnergyBalanceModel.

import StyledStrings, Statistics

export Progress, update!
export crossmean, condcopy!, condcopy, zeroref!

# progress bar
mutable struct Progress
    title::String
    total::Int
    current::Int
    last::Int # last printed progress
    started::Float64 # start time
    updated::Float64 # last external update time
    freq::Float64 # external update frequency
    infofeed::Function # Function(done::Bool, args...)::String
    width::Int # number of characters wide, including progress texts
    barwidth::Int # width of the progress bar
    lines::Int # lines printed
    runners::Tuple{Vararg{Char,4}} # characters to use as runners
    updates::Int # number of external updates

    function Progress(
        total::Int,
        title::String="Progress", freq::Float64=1.0, infofeed::Function=(done::Bool -> "");
        width::Int=50
    )
        barwidth = width - (ndigits(total) * 2 + 1) - 2 - 5 - 3 # current/total [=> ] xx.x%
        return new(
            title,
            total,
            -1, # current
            0, # last
            NaN, # started
            NaN, # updated
            freq,
            infofeed,
            width,
            barwidth,
            0, # lines
            ('◓', '◑', '◒', '◐'), # runners
            0 # updates
        ) # new
    end # function Progress
end # struct Progress

function display_time(time::Float64)::String
    if isfinite(time) # remaining time unknown
        timeint = round(Int, time)
        min = fld(timeint, 60)
        sec = timeint % 60
        str = string(min, ':', string(sec, pad=2))
    else # !isfinite(time)
        str = "-:--"
    end # if isfinite, else
    return str
end # function display_time

function output!(prog::Progress, feedargs::Tuple{Vararg{Any}}=())::Nothing
    now = time()
    isdone = false
    # avoid over-updating
    if prog.current > prog.total
        return nothing
    end # if >
    # clear previous lines
    while prog.lines > 0
        print("\033[A\033[2K") # move up one line and clear the line
        prog.lines -= 1 # !
    end # while >
    # title
    println(StyledStrings.styled"{bold,region,warning:$(prog.title)}")
    prog.lines += 1 # !
    # get bar and info strings
    elapsed = display_time(now - prog.started)
    if prog.current >= prog.total # done
        isdone = true
        # progress
        barstr = StyledStrings.annotatedstring(
            # current/total
            lpad(StyledStrings.styled"{success:$(prog.current)}", ndigits(prog.total) + 1), '/', prog.total,
            # bar
            " [", StyledStrings.styled"""{success:$(repeat("=", prog.barwidth))}""", "] ",
            # percentage
            lpad(StyledStrings.styled"{success:$(round(Int, prog.current/prog.total*100))%}", 5)
        ) # StyledStrings.annotatedstring
        # time and speed
        speed = prog.current / (now - prog.started)
        togo = display_time((prog.total - prog.current) / speed)
        prompt = StyledStrings.styled"{success:{bold:Done} ✓}"
    else # in progress
        # progress
        done = floor(Int, prog.current / prog.total * prog.barwidth) # number of chars to fill =
        barstr = StyledStrings.annotatedstring(
            # current/total
            lpad(StyledStrings.styled"{info:$(prog.current)}", ndigits(prog.total) + 1), '/', prog.total,
            # bar
            " [",
            StyledStrings.styled"""{info:$(repeat("=", done))>}""",
            repeat(" ", max(prog.barwidth - done - 1, 0)),
            "] ",
            # percentage
            lpad(StyledStrings.styled"{info:$(round(prog.current/prog.total*100, digits=1))%}", 5)
        ) # StyledStrings.annotatedstring
        # time and speed
        speed = (prog.current - prog.last) / (now - prog.updated)
        togo = display_time((prog.total - prog.current) / speed)
        prompt = StyledStrings.styled"{info:{bold:In progress} $(prog.runners[mod1(prog.updates, 4)])}"
    end # if >=, else
    prog.last = prog.current # !
    prog.updated = now # !
    prog.updates += 1 # !
    if !isfinite(speed) # no speed info
        spdstr = "-/sec"
    elseif (speed >= 1.0) || (iszero(speed)) # speed > 1.0
        spdstr = string(round(speed, digits=2), "/sec")
    else # speed < 1.0
        spdstr = string(round(1.0 / speed, digits=2), "sec/1")
    end # if >, elseif, else
    timespeed = StyledStrings.annotatedstring(
        ' ',
        StyledStrings.styled"{$(isdone ? :success : :info):$elapsed}",
        "/",
        StyledStrings.styled"{note:-$togo}",
        ' ',
        spdstr
    ) # StyledStrings.annotatedstring
    infopaddings = repeat(" ", max(prog.width - length(timespeed) - length(prompt), 1))
    # output bar and info
    println(barstr)
    prog.lines += 1 # !
    println(timespeed, infopaddings, prompt)
    prog.lines += 1 # !
    # update user custom info
    userstr::String = prog.infofeed(isdone, feedargs...)
    userstrvec = split(userstr)
    annotatedvec = map((s::String -> StyledStrings.styled" {note:$s}"), userstrvec)
    foreach(s::Base.AnnotatedString{String} -> println(s), annotatedvec)
    prog.lines += length(annotatedvec) # !
    return nothing
end # function output

function update!(prog::Progress, current::Int=prog.current+1, feedargs::Tuple{Vararg{Any}}=())::Nothing
    # internal update
    prog.current = current # !
    # initialise if not started
    if isnan(prog.started)
        prog.started = time() # !
        prog.updated = time() - prog.freq # force immediate external update # !
    end # if isnan
    # external update
    if (prog.current >= prog.total) || (time() - prog.updated >= prog.freq)
        output!(prog, feedargs)
    end # if ||
    return nothing
end # function update!

@inline function crossmean(vecvec::Vector{Vector{T}})::Vector{T} where T<:Number
    @boundscheck if !all(length.(vecvec) .== length(vecvec[1]))
        throw(BoundsError("All vectors must be the same length."))
    end # if !
    return map((xi::Int -> Statistics.mean([vecvec[ti][xi] for ti in eachindex(vecvec)])), eachindex(vecvec[1]))
end # function crossmean

# conditional copy in place
@inline function condcopy!(to::Vector{T}, from::T, cond::Function, ref::Vector{T}=to)::Vector{T} where T
    @. to[cond(ref)] = from # !
    return to
end # function condcopy!

@inline (condcopy(to::Vector{T}, from::T, cond::Function, ref::Vector{T}=to)::Vector{T}) where T =
    condcopy!(copy(to), from, cond, ref)

# replace entries with zeros in ref with zeros in place in v
@inline (zeroref!(v::Vector{T}, ref::Vector{T})::Vector{T}) where T<:Number = condcopy!(v, zero(T), iszero, ref)

end # module Utilities


module Infrastructure # EnergyBalanceModel.

using ..Utilities

import EnergyBalanceModel
import SparseArrays, Statistics

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
(Base.keys(coll::Collection{V})::AbstractSet{Symbol}) where V = keys(getfield(coll, :dict))

struct SpaceTime{F<:Function}
    nx::Int # number of evenly spaced latitudinal gridboxes (equator to pole)
    x::Vec # grid
    xmodifier::F # function to modify x grid
    dur::Int # duration of simulation in years
    nt::Int # number of timesteps per year (limited by numerical stability)
    dt::Float64 # timestep
    t::Vec # time vector in a year
    T::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64} # full time series
    winter::@NamedTuple{t::Float64, inx::Int}
    summer::@NamedTuple{t::Float64, inx::Int}

    function SpaceTime{F}(
        xrange::Tuple{Float64,Float64}, nx::Int, nt::Int, dur::Int;
        winter::Float64=0.26125, summer::Float64=0.77375
    ) where F<:Function
        dx = (xrange[2]-xrange[1]) / nx
        x = F.instance.(collect(dx/2.0 : dx : xrange[2] - dx/2.0))
        xmodifier = F.instance
        dt = 1.0 / nt
        t = collect(range(dt/2.0, 1.0 - dt/2.0, nt))
        T = dt/2.0 : dt : dur - dt/2.0
        winterinx = round(Int, nt*winter)
        summerinx = round(Int, nt*summer)
        return new{F}(
            nx, x, xmodifier, dur, nt, dt, t, T, (t=winter, inx=winterinx), (t=summer, inx=summerinx)
        )
    end # function SpaceTime
end # struct SpaceTime{F<:Function}

SpaceTime(
    ::typeof(identity), nx::Int, nt::Int, dur::Int; winter::Float64=0.26125, summer::Float64=0.77375
) = SpaceTime{typeof(identity)}((0.0, 1.0), nx, nt, dur, winter=winter, summer=summer)
SpaceTime(
    ::typeof(sin), nx::Int, nt::Int, dur::Int; winter::Float64=0.26125, summer::Float64=0.77375
) = SpaceTime{typeof(sin)}((0.0, pi/2.0), nx, nt, dur, winter=winter, summer=summer)
SpaceTime(
    nx::Int, nt::Int, dur::Int; winter::Float64=0.26125, summer::Float64=0.77375
) = SpaceTime(identity, nx, nt, dur, winter=winter, summer=summer)

struct Forcing
    flat::Bool # forcing is always at base
    base::Float64 # base forcing
    peak::Float64 # peak forcing
    cool::Float64 # forcing after cooldown
    holdyrs::Tuple{Int,Int} # years to hold at (base, peak) forcing
    rates::Tuple{Float64,Float64} # rates of change
    domain::Tuple{Vararg{Int,6}} # years at which forcing pattern changes

    # constant forcing
    Forcing(base::Float64) = new(
        true, base, base, base, (0, 0), (0.0, 0.0), (0, 0, 0, 0, 0, 0)
    )
    # warming/cooling forcing
    function Forcing(
        base::Float64, peak::Float64, cool::Float64, holdyrs::Tuple{Int,Int}, rates::Tuple{Float64,Float64}
    )
        domainvec = zeros(Int, 6)
        # hold at base
        @. domainvec[2:6] += holdyrs[1]
        # time to warm
        warming = (peak - base) / rates[1]
        isinteger(warming) ?
            @.(domainvec[3:6] += warming) :
            throw(ArgumentError("Warming time must be integer. Got $warming y."))
        # hold at peak
        @. domainvec[4:6] += holdyrs[2]
        # time to cool
        cooling = (peak - cool) / rates[2]
        isinteger(cooling) ?
            @.(domainvec[5:6] += cooling) :
            throw(ArgumentError("Cooling time must be integer. Got $cooling y."))
        # hold at cool
        domainvec[6] = Inf
        return new(false, base, peak, cool, holdyrs, rates, Tuple(domainvec))
    end # function Forcing
end # struct Forcing

# evaluate forcing at time T (in years)
function (forcing::Forcing)(T::Float64)::Float64
    if forcing.flat
        f = forcing.base
    else # !forcing.flat
        if T < forcing.domain[2] # hold at base
            f = forcing.base
        elseif T < forcing.domain[3] # warming
            f = forcing.base + forcing.rates[1] * (T-forcing.domain[2])
        elseif T < forcing.domain[4] # hold at peak
            f = forcing.peak
        elseif T < forcing.domain[5] # cooling
            f = forcing.peak - forcing.rates[2] * (T-forcing.domain[4])
        else # hold at cool
            f = forcing.cool
        end # if <, elseif*3, else
    end # if forcing.flat, else
    return f
end # function (forcing::Forcing)

struct Solutions
    spacetime::SpaceTime{<:Function} # space and time which solutions are defined on
    ts::Vec # time vector for stored solution
    forcing::Forcing # climate forcing
    parameters::Collection{Float64} # model parameters
    initconds::Collection{Vec} # initial conditions
    lastonly::Bool # store only last year of solution
    debug::Expr # store debug variables
    raw::Collection{Vector{Vec}} # solution storage
    seasonal::@NamedTuple{
        winter::Collection{Vector{Vec}}, summer::Collection{Vector{Vec}}, avg::Collection{Vector{Vec}}
    } # seasonal peak and annual avg

    function Solutions(
        st::SpaceTime{<:Function}, forcing::Forcing, par::Collection{Float64}, init::Collection{Vec}, vars::Set{Symbol},
        lastonly::Bool=true;
        debug::Expr=Expr(:block)
    ) # Solutions
        if lastonly
            dur_store = 1
            ts::Vec = st.dur-1.0 + st.dt/2.0 : st.dt : st.dur - st.dt/2.0
        else # !lastonly
            dur_store = st.dur
            ts = st.dt/2.0 : st.dt : st.dur - st.dt/2.0
        end # if lastonly, else
        if debug != Expr(:block)
            push!(vars, :debug)
        end # if !=
        # construct raw solution storage
        solraw = Collection{Vector{Vec}}()
        foreach((var::Symbol -> setproperty!(solraw, var, Vector{Vec}(undef, length(ts)))), vars)
        # construct seasonal solution storage template
        seasonaltemp = Collection{Vector{Vec}}()
        foreach((var::Symbol -> setproperty!(seasonaltemp, var, Vector{Vec}(undef, st.dur))), vars)
        return new(
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
end # struct Solutions

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
    model == :MIZ ? default_parameters(miz_paramset) : default_parameters(classic_paramset)

# calculate diffusion operator matrix
let diffop::SparseArrays.SparseMatrixCSC{Float64,Int64} = SparseArrays.spzeros(Float64, 0, 0)
    @inline function get_diffop(nx::Int)::SparseArrays.SparseMatrixCSC{Float64,Int64}
        if size(diffop) !== (nx, nx) # recalculate diffusion operator
            dx = 1.0 / nx
            xb = dx : dx : 1.0-dx
            lambda = @. (1 - xb^2) / dx^2
            l1 = pushfirst!(-copy(lambda), 0.0)
            l2 = push!(-copy(lambda), 0.0)
            l3 = -l1 - l2
            diffop = SparseArrays.spdiagm(-1 => -l1[2:nx], 0 => -l3, 1 => -l2[1:nx-1])
        end
        return diffop
    end # function get_diffop

    @eval (@__MODULE__) get_diffop(nx::Int)::SparseArrays.SparseMatrixCSC{Float64,Int64} = $get_diffop(nx)
end # let diffop

# diffusion for equal spaced grid
@inline (diffusion!(
    base::VT, T::VT, st::SpaceTime{typeof(identity)}, par::Collection{Float64}
)::VT) where VT<:Vector{<:Number} = base .+= par.D * get_diffop(st.nx) * T

# diffusion for non-equal spaced grid
@inline function diffusion!(
    base::VT, T::VT, st::SpaceTime{<:Function}, par::Collection{Float64}
)::VT where VT<:Vector{<:Number}
    diffT = diff(T)
    diffx = diff(st.x)
    i = 2 : st.nx-1
    xi = @view st.x[i]
    xim = @view st.x[i.-1]
    xip = @view st.x[i.+1]
    xxph = @. (xip + xi) / 2.0
    xxmh = @. (xi + xim) / 2.0
    @inbounds @. base[i] +=
        par.D * ((1-xxph^2) * diffT[i]/diffx[i] - (1-xxmh^2) * diffT[i-1]/diffx[i-1]) /
        (xxph - xxmh) # !
    # base[1] += par.D * (-2.0)*st.x[1] * diffT[1] / diffx[1] # !
    base[1] += par.D * diffT[1] / diffx[1]^2 # ! # TODO why?
    base[end] += par.D * (-2.0)*st.x[end] * diffT[end] / diffx[end] # !
    return base
end

@inline diffusion(T::Vec, st::SpaceTime{<:Function}, par::Collection{Float64})::Vec =
    diffusion!(zeros(Float64, length(T)), T, st, par)

const D∇² = diffusion
const D∇²! = diffusion!

# calculate annual mean
function annual_mean(annusol::Solutions)::Collection{Vec}
    # calculate annual mean for each variable except temperatures
    means = Collection{Vec}()
    foreach(
        (var::Symbol -> setproperty!(means, var, crossmean(getproperty(annusol.raw, var)))),
        keys(annusol.raw)
    )
    return means
end # function annual_mean

annual_mean(forcing::Forcing, st::SpaceTime{<:Function}, year::Int)::Float64 = Statistics.mean(forcing.(year-1 .+ st.t))

function savesol!(
    sols::Solutions, annusol::Solutions, vars::Collection{Vec}, tinx::Int
)::Solutions
    varscp = deepcopy(vars) # avoid reference issues
    year = ceil(Int, sols.spacetime.T[tinx])
    ti = mod1(tinx, sols.spacetime.nt) # index of time in the year
    # save raw data to annual
    foreach(keys(annusol.raw)) do var::Symbol
        getproperty(annusol.raw, var)[ti] = getproperty(varscp, var) # !
    end # foreach do
    # save raw data
    if !sols.lastonly # save all raw data
        foreach(
            (var::Symbol -> setindex!(getproperty(sols.raw, var), getproperty(varscp, var), tinx)),
            keys(sols.raw)
        )
    elseif tinx > length(sols.spacetime.T) - sols.spacetime.nt # save the raw data of the last year
        foreach(
            (var::Symbol -> setindex!(getproperty(sols.raw, var), getproperty(varscp, var), ti)),
            keys(sols.raw)
        )
    end # if !, elseif
    # save seasonal data
    if ti == sols.spacetime.winter.inx
        foreach(
            (var::Symbol -> setindex!(getproperty(sols.seasonal.winter, var), getproperty(varscp, var), year)),
            keys(sols.seasonal.winter)
        )
    elseif ti == sols.spacetime.summer.inx
        foreach(
            (var::Symbol -> setindex!(getproperty(sols.seasonal.summer, var), getproperty(varscp, var), year)),
            keys(sols.seasonal.summer)
        )
    elseif ti == sols.spacetime.nt # calculate annual average
        means = annual_mean(annusol)
        foreach(
            (var::Symbol -> setindex!(getproperty(sols.seasonal.avg, var), getproperty(means, var), year)),
            keys(sols.seasonal.avg)
        )
    end # if ==, elseif*2
    return sols
end # function savesol!

function integrate(
    model::Symbol, st::SpaceTime{<:Function}, forcing::Forcing, par::Collection{Float64}, init::Collection{Vec};
    lastonly::Bool=true, debug::Expr=Expr(:block)
)::Solutions
    # initialise
    vars = deepcopy(init)
    solvars = Set{Symbol}((:E, :T)) # always solve for these
    if model === :MIZ # add MIZ variables
        union!(solvars, Set{Symbol}((:Ei, :Ew, :Ti, :Tw, :h, :D, :phi, :n)))
    end
    Modu::Module = EnergyBalanceModel.eval(model)
    sols = Solutions(st, forcing, par, init, solvars, lastonly, debug=debug)
    annusol = Solutions(st, forcing, par, init, solvars, true, debug=debug) # for calculating annual means
    progress = Progress(length(st.T), "Integrating")
    update!(progress)
    # loop over time
    for ti in eachindex(st.T)
        Modu.step!(st.t[mod1(ti, st.nt)], forcing(st.T[ti]), vars, st, par, debug=debug)
        savesol!(sols, annusol, vars, ti)
        update!(progress)
    end # for in
    return sols
end # function integrate

end # module Infrastructure


module MIZ # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import NonlinearSolve

# solar radiation absorbed on ice and water
@inline (solar!(base::VT, x::Vec, t::Float64, ::Val{true}, par::Collection{Float64})::VT) where VT<:Vector{<:Number} =
    @. (base += par.ai * (par.S0 - par.S1 * x * cos(2.0*pi * t) - par.S2 * x^2))
@inline (solar!(base::VT, x::Vec, t::Float64, ::Val{false}, par::Collection{Float64})::VT) where VT<:Vector{<:Number} =
    @. (base += (par.a0 - par.a2 * x^2) * (par.S0 - par.S1 * x * cos(2.0*pi * t) - par.S2 * x^2))

@inline solar(x::Vec, t::Float64, ice::Bool, par::Collection{Float64})::Vec = solar!(
    zeros(Float64, length(x)), x, t, Val(ice), par
)

# temperatures
@inline function Tbar!(Ti::VT, Tw::Vec, phi::Vec)::VT where VT<:Vector{<:Number}
    Ti .*= phi # !
    @. Ti += (1 - phi) * Tw # !
    return Ti
end
@inline Tbar(Ti::Vec, Tw::Vec, phi::Vec)::Vec = Tbar!(copy(Ti), Tw, phi)
const T̄ = Tbar
const T̄! = Tbar!

@inline water_temp(Ew::Vec, phi::Vec, par::Collection{Float64})::Vec = @. par.Tm + Ew / ((1-phi) * par.cw)
@inline (ice_temp(T0::VT, par::Collection{Float64})::VT) where VT<:Vector{<:Number} = min.(T0, par.Tm)

function T0eq(
    T0::VT,
    args::@NamedTuple{
        x::Vec, t::Float64, h::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime{F}, par::Collection{Float64}
    }
)::VT where {VT<:Vector{<:Number}, F<:Function} # T0eq
    vec = @. args.par.k * (args.par.Tm - T0) / args.h # SCM
    solar!(vec, args.x, args.t, Val(true), args.par) # solar on ice
    @. vec += (-args.par.A) - args.par.B * (T0 - args.par.Tm) # OLR
    diffusion!(vec, Tbar!(ice_temp(T0, args.par), args.Tw, args.phi), args.st, args.par) # diffusion
    vec .+= args.f # forcing
    return vec
end # function T0eq

let T0::Vec = zeros(Float64, 100) # let T0 be a persistent variable
    function solveTi(
        x::Vec, t::Float64, h::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime{<:Function}, par::Collection{Float64}
    )::Vec
        h1 = condcopy(h, 1.0, iszero) # avoid division by zero when solving T0
        if length(T0) != length(x)
            T0 = zeros(Float64, length(x)) # initialise T0
        end # if !=
        T0sol = NonlinearSolve.solve(
            NonlinearSolve.NonlinearProblem(T0eq, T0, (; x, t, h=h1, Tw, phi, f, st, par)),
            NonlinearSolve.TrustRegion()
        )
        # TODO test for solving failure
        T0 = T0sol.u
        Ti = ice_temp(T0, par)
        zeroref!(Ti, h) # set Ti to 0 where no ice
        return Ti
    end # function solveTi

    @eval (@__MODULE__) solveTi(
        x::Vec, t::Float64, h::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime{<:Function}, par::Collection{Float64}
    )::Vec = $solveTi(x, t, h, Tw, phi, f, st, par)
end # let T0

# lateral melt rate
wlat(Tw::Vec, par::Collection{Float64})::Vec = @. par.m1 * (Tw - par.Tm^par.m2)

# concentration
function concentration(Ei::Vec, h::Vec, par::Collection{Float64})::Vec
    phi = @. -Ei / (par.Lf * h)
    zeroref!(phi, h)
    condcopy!(phi, 1.0, >(1.0)) # correct concentration
    # phi = @. Float64(Ei.<0.0) # reproducing WE15
    return phi
end # function concentration

# floe number
function num(D::Vec, phi::Vec, par::Collection{Float64})::Vec
    n = @. phi / (par.alpha * D^2)
    zeroref!(n, D)
    return n
end # function num

# lead region area
function area_lead(D::Vec, phi::Vec, n::Vec, par::Collection{Float64})::Vec
    ring = @. par.alpha * n * ((D + 2.0*par.rl)^2 - D^2)
    return min.(ring, 1.0.-phi)
end # function area_lead

# fluxes
function vert_flux(
    x::Vec, t::Float64, ice::Bool, Ti::Vec, Tw::Vec, phi::Vec, f::Float64,
    st::SpaceTime{<:Function}, par::Collection{Float64}
)::Vec
    L = @. par.A + par.B * ($Tbar(Ti, Tw, phi) - par.Tm) # OLR
    return @. $solar(x, t, ice, par) - L + $(diffusion(Tbar(Ti, Tw, phi), st, par)) + par.Fb + f
end # function vert_flux

function lat_flux(h::Vec, D::Vec, Tw::Vec, phi::Vec, par::Collection{Float64})::Vec
    Flat = @. phi * h * par.Lf * $wlat(Tw, par) * pi / (par.alpha * D)
    zeroref!(Flat, D)
    return Flat
end # function lat_flux

function redistributeE(rEi::Vec, rEw::Vec)::@NamedTuple{Ei::Vec, Ew::Vec, psiEidt::Vec, psiEwdt::Vec}
    cEi = clamp.(rEi, -Inf, 0.0)
    cEw = clamp.(rEw, 0.0, Inf)
    psiEidt = rEi .- cEi # +
    psiEwdt = rEw .- cEw # -
    Ei = cEi .+ psiEwdt # -
    Ew = cEw .+ psiEidt # +
    return (; Ei, Ew, psiEidt, psiEwdt)
end # function redistributeE

# redistribution functions
function split_psiEw(psiEw::Vec, phi::Vec, Al::Vec)::@NamedTuple{Ql::Vec, Qp::Vec}
    Ql = @. Al / (1-phi) * psiEw
    condcopy!(Ql, 0.0, isone, phi) # fix rounding errors
    Qp = psiEw - Ql
    return (; Ql, Qp)
end # function split_psiEw

psinplus(Qp::Vec, par::Collection{Float64})::Vec = -Qp / (par.Lf * par.alpha * par.Dmin^2 * par.hmin)

function average(f::Vec, fn::Float64, n::Vec, dn::Vec)::Vec
    total = n .+ dn
    avgd = @. (n*f + dn*fn) / total
    zeroref!(avgd, total)
    return avgd
end # function average

# differential equations
Ei_t(phi::Vec, Fvi::Vec, Flat::Vec)::Vec = @. phi * Fvi + Flat
Ew_t(phi::Vec, Fvw::Vec, Flat::Vec)::Vec = @. (1 - phi) * Fvw - Flat
h_t(Fvi::Vec, par::Collection{Float64})::Vec = -1/par.Lf * Fvi
function D_t(h::Vec, D::Vec, Tw::Vec, phi::Vec, Ql::Vec, par::Collection{Float64})::Vec
    lat_melt = -pi / 2.0*par.alpha * wlat(Tw, par)
    lat_grow = @. -D / (2 * par.Lf * h * phi) * Ql
    weld = @. par.kappa * par.alpha / 4 * phi * D^3
    zeroref!(lat_grow, h)
    return @. lat_melt + lat_grow + weld
end # function D_t

forward_euler(var::Vec, grad::Vec, dt::Float64)::Vec = @. var + grad * dt

function step!(
    t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime{<:Function}, par::Collection{Float64};
    debug::Expr=Expr(:block)
)::Collection{Vec}
    # update temperature
    vars.Tw = water_temp(vars.Ew, vars.phi, par) # !
    condcopy!(vars.Tw, 0.0, isnan) # eliminate NaNs for calculations
    vars.Ti = solveTi(st.x, t, vars.h, vars.Tw, vars.phi, f, st, par) # !
    # update floe number
    vars.n = num(vars.D, vars.phi, par) # !
    # calculate fluxes
    Fvi = vert_flux(st.x, t, true, vars.Ti, vars.Tw, vars.phi, f, st, par)
    Fvw = vert_flux(st.x, t, false, vars.Ti, vars.Tw, vars.phi, f, st, par)
    Flat = lat_flux(vars.h, vars.D, vars.Tw, vars.phi, par)
    # update enthalpy
    rEi = forward_euler(vars.Ei, Ei_t(vars.phi, Fvi, Flat), st.dt)
    rEw = forward_euler(vars.Ew, Ew_t(vars.phi, Fvw, Flat), st.dt)
    Epsidt = redistributeE(rEi, rEw)
    vars.Ei = Epsidt.Ei # !
    vars.Ew = Epsidt.Ew # !
    # update floe size and thickness
    Al = area_lead(vars.D, vars.phi, vars.n, par)
    Qlp = split_psiEw(Epsidt.psiEwdt/st.dt, vars.phi, Al)
    dn = st.dt * psinplus(Qlp.Qp, par) # number of new pancakes
    rD = forward_euler(vars.D, D_t(vars.h, vars.D, vars.Tw, vars.phi, Qlp.Ql, par), st.dt)
    vars.D = average(rD, par.Dmin, vars.n, dn) # new pancakes # !
    clamp!(vars.D, par.Dmin, par.Dmax)
    zeroref!(vars.D, vars.Ei) # correct round off errors
    rh = forward_euler(vars.h, h_t(Fvi, par), st.dt)
    clamp!(rh, 0.0, Inf) # avoid overshooting to negative thickness
    vars.h = average(rh, par.hmin, vars.n, dn) # new pancakes # !
    # update concentration
    vars.phi = concentration(vars.Ei, vars.h, par) # !
    # update total energy and temperature
    zeroref!(vars.Ei, vars.h) # correct round off errors
    vars.E = @. vars.phi * vars.Ei + (1 - vars.phi) * vars.Ew # !
    vars.T = Tbar(vars.Ti, vars.Tw, vars.phi) # !
    # debug
    if debug != Expr(:block)
        vars.debug = eval(debug) # !
    end # if !=
    # set NaNs to no existence
    condcopy!(vars.Ti, NaN, iszero, vars.Ei)
    condcopy!(vars.Tw, NaN, isone, vars.Ew)
    return vars
end # function step

end # module MIZEBM


module Classic # EnergyBalanceModel.

using ..Infrastructure

using AppleAccelerate
import LinearAlgebra, SparseArrays

let id::UInt64 = UInt64(0),
    cg_tau::Float64 = NaN,
    dt_tau::Float64 = NaN,
    dc::Float64 = NaN,
    kappa::Matrix{Float64} = Matrix{Float64}(undef, 100, 100),
    S::Matrix{Float64} = Matrix{Float64}(undef, 100, 2001),
    M::Float64 = NaN,
    aw::Vec = Vec(undef, 100),
    kLf::Float64 = NaN # let ,*8

    @inline function get_statics(st::SpaceTime{<:Function}, par::Collection{Float64})::@NamedTuple{
        cg_tau::Float64, dt_tau::Float64, dc::Float64, kappa::Matrix{Float64},
        S::Matrix{Float64}, M::Float64, aw::Vec, kLf::Float64
    }
        if id != hash((st, par)) # recompute only if st or par changed
            # Difinitions for implicit scheme for Tg
            cg_tau = par.cg / par.tau
            dt_tau = st.dt / par.tau
            dc = dt_tau * cg_tau
            kappa = (1+dt_tau) * LinearAlgebra.I(st.nx) - st.dt * par.D * get_diffop(st.nx) / par.cg
            # Seasonal forcing [WE15 Eq. (3)]
            S = repeat(par.S0 .- par.S2 * st.x.^2, 1, st.nt) -
                repeat(par.S1 * cos.(2.0*pi*st.t'), st.nx, 1) .* repeat(st.x, 1, st.nt)
            S = hcat(S, S[:,1])
            # Further definitions
            M = par.B + cg_tau
            aw = @. par.a0 - par.a2 * st.x^2
            kLf = par.k * par.Lf
            # update id
            id = hash((st, par))
        end # if !=
        return (; cg_tau, dt_tau, dc, kappa, S, M, aw, kLf)
    end # function get_statics

    @eval (@__MODULE__) get_statics(st::SpaceTime{<:Function}, par::Collection{Float64})::@NamedTuple{
        cg_tau::Float64, dt_tau::Float64, dc::Float64, kappa::Matrix{Float64},
        S::Matrix{Float64}, M::Float64, aw::Vec, kLf::Float64
    } = $get_statics(st, par) # @eval
end # let id, cg_tau, dt_tau, dc, kappa, S, M, aw, kLf

function step!(
    t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime{<:Function}, par::Collection{Float64};
    debug::Expr=Expr(:block)
)::Collection{Vec}
    # get static variables
    stat = get_statics(st, par)
    # get time index
    i = round(Int, mod1((t + st.dt/2.0) * st.nt, st.nt))
    # forcing
    alpha = @. stat.aw * (vars.E>0) + par.ai * (vars.E<0) # WE15 Eq. (4)
    C = @. alpha*stat.S[:,i] + stat.cg_tau*vars.Tg - par.A + f
    # surface temperature
    T0 = @. C / (stat.M - stat.kLf/vars.E) # WE15 Eq. (A3)
    vars.T = @. vars.E/par.cw * (vars.E>=0) + T0 * (vars.E<0)*(T0<0) # WE15 Eq. (9)
    # Forward Euler for E
    @. vars.E += st.dt * (C - stat.M*vars.T + par.Fb) # WE15 Eq. (A2)
    # Implicit Euler for Tg
    vars.Tg =
        (stat.kappa - SparseArrays.spdiagm(stat.dc ./ (stat.M .- stat.kLf./vars.E) .* (T0.<0).*(vars.E.<0))) \
        (
            vars.Tg +
            (
                stat.dt_tau * (vars.E/par.cw.*(vars.E.>=0) +
                (par.ai*view(stat.S, :, i+1) .- par.A .+ f) ./ (stat.M .- stat.kLf./vars.E) .* (T0.<0).*(vars.E.<0))
            )
        ) # () # vars.Tg # WE15 Eq. (A1)
    # debug
    if debug != Expr(:block)
        vars.debug = eval(debug) # !
    end # if !=
    return vars
end # function step

end # module ClassicEBM


module Plot # EnergyBalanceModel.

end # module Plot


using .Infrastructure

export Vec, Collection, SpaceTime, Forcing, Solutions
export miz_paramset, classic_paramset, default_parameters
export integrate

end # module EnergyBalanceModel
