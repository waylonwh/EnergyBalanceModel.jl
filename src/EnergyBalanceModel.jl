"""
    EnergyBalanceModel

A comprehensive package for solving a classic Energy Balance Model (EBM) (Wagner and
Eisenman, 2015) and an extended EBM with the inclusion of a Marginal Ice Zone (MIZ). Other
utilities for data handling and visualization are also provided.

To get started, define a space-time domain, a forcing function, parameters, and initial
conditions. Then call `integrate` to run the model. See documentation of `SpaceTime`,
`Forcing`, `default_parameters`, `Collection`, and `integrate` for details. The following
example runs the EBM with MIZ for 30 years on a 180-point latitudinal grid equally spaced in
sine latitude, with 2000 timesteps per year and a constant forcing of 0.0. Then the results
are saved and plotted.

```julia-repl
julia> using EnergyBalanceModel

julia> st = SpaceTime{sin}(180, 2000, 30)
SpaceTime{sin} with:
  180 latitudinal gridboxes: [0.00436331, 0.0130896, … 762, 0.999914, 0.99999]
  2000 timesteps per year: [0.00025, 0.00075, 0.001 … 99875, 0.99925, 0.99975]
  30 years of simulation: t∈[0,30]
  winter at t=0.26125, summer at t=0.77375

julia> forcing = Forcing(0.0)
Forcing{true}(0.0) is constant:
  F(t)=0.0, t∈[0,∞)

julia> par = default_parameters(:MIZ)
Collection{Float64} with 22 entries:
  :Dmax  => 156.0
  :a2    => 0.1
  :alpha => 0.66
  :m1    => 50.4576
  :D     => 0.6
  :S1    => 338.0
  :B     => 2.1
  :cw    => 9.8
  :rl    => 0.5
  :Fb    => 4.0
  ⋮      => ⋮

julia> init = Collection{Vec}(
           :Ei => zeros(st.nx),
           :Ew => zeros(st.nx),
           :h => zeros(st.nx),
           :D => zeros(st.nx),
           :phi => zeros(st.nx)
       )
Collection{Vector{Float64}} with 5 entries:
  :Ei  => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, …
  :D   => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, …
  :h   => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, …
  :phi => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, …
  :Ew  => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, …

julia> sols = integrate(:MIZ, st, forcing, par, init)
Integrating
 60000/60000 [━━━━━━━━━━━━━━━━━━━━━━━━━━━━━]  100%
 1:57/-0:00 511.24/sec                      Done ✓
 t = 30.0
Solutions{sin, true} with:
  10 solution variables: [:T, :Ei, :Ti, :D, :n, :h, :phi, :E, :Ew, :Tw]
  on 180 latitudinal gridboxes: [0.00436331, 0.0130896 … 2, 0.999914, 0.99999]
  and 2000 timesteps: 29.00025:0.0005:29.99975
  with forcing Forcing{true}(0.0) (constant forcing)

julia> save(sols, "./miz_sol.jld2")
"./miz_sol.jld2"

julia> plot_raw(sols)
```

See the documentation for submodules `IO` and `Plot` for details on data handling and
visualization.
"""
module EnergyBalanceModel # .


module Utilities # EnergyBalanceModel.

import UUIDs, StyledStrings as SS, Statistics as Stats, TimeZones as TZ

export Progress, update!
export safehouse, house!, favorite!, note!, retrieve
export @persistent, iobuffer, unique_id, reprhex
export crossmean, hemispheric_mean
export condset!, condset, zeroref!

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
    runners::NTuple{4,Char} # characters to use as runners
    updates::Int # number of external updates

    function Progress(
        total::Int,
        title::String="Progress", freq::Float64=1.0;
        width::Int=50, infofeed::Function=(_ -> "")
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

# members in the safehouse
mutable struct Refugee{M,T}
    varname::Symbol
    id::UInt32
    housed::TZ.ZonedDateTime
    val::T
    note::String

    function Refugee{M}(var::Symbol) where M
        val = deepcopy(getproperty(M, var))
        return new{M,typeof(val)}(var, unique_id(), TZ.now(TZ.localzone()), val, "")
    end # function Refugee{M}
end # struct Refugee{M,T}

function Base.show(io::IO, refugee::Refugee{M,T})::Nothing where {M, T}
    print(
        io,
        typeof(refugee), '(', refugee.varname, '#', reprhex(refugee.id), " = "
    )
    show(io, refugee.val)
    print(io, ')')
end # function Base.show

function Base.show(io::IO, ::MIME"text/plain", refugee::Refugee{M,T})::Nothing where {M, T}
    println(
        io,
        typeof(refugee), '(', refugee.varname, '#', reprhex(refugee.id), ')', " housed at ", refugee.housed, ':'
    )
    buffer = iobuffer(io; sizemodifier=(0, -2))
    show(buffer, MIME("text/plain"), refugee.val)
    str = String(take!(buffer.io))
    print(io, string("  ", replace(str, '\n' => "\n  ")))
    return nothing
end # function Base.show

# safehouse to hold results before being overwritten
struct Safehouse{M}
    variables::Dict{Symbol,Vector{UInt32}}
    favorites::Dict{Symbol,UInt32}
    refugees::Dict{UInt32,Refugee{M}}

    function Safehouse{M}(name::Symbol=:SAFEHOUSE) where M
        safehouse = new{M}(Dict{Symbol,Vector{UInt32}}(), Dict{Symbol,UInt32}(), Dict{UInt32,Refugee{M}}())
        @eval M const $name = $safehouse
        return safehouse
    end # function Safehouse{M}
end # struct Safehouse{M}

(Base.show(io::IO, safehouse::Safehouse{M})::Nothing) where M = print(
    io,
    typeof(safehouse),
    '(',
    join([string(length(safehouse.variables[v]), '@', v) for v in keys(safehouse.variables)], ", "),
    ')'
)

function Base.show(io::IO, ::MIME"text/plain", safehouse::Safehouse{M})::Nothing where M
    print(
        io,
        typeof(safehouse), " with ", length(safehouse.refugees), " refugees in ",
        length(safehouse.variables), " variables:"
    )
    for ids in values(safehouse.variables), id in ids
        print(io, "\n  ")
        show(io, safehouse.refugees[id])
    end # for ids, id
    return nothing
end # function Base.show

macro persistent(exprs...)
    # syntax tree operations
    findexpr(_, ::Symbol)::Nothing = nothing
    function findexpr(expr::Expr, head::Symbol)::Union{Expr,Nothing}
        if expr.head === head
            return expr
        elseif isempty(expr.args)
            return nothing
        else
            for arg in expr.args
                funcexpr = findexpr(arg, head)
                if !isnothing(funcexpr)
                    return funcexpr
                end
            end
            return nothing
        end # if ==
    end # function findexpr
    sign2call(expr::Symbol)::Symbol = expr
    function sign2call(expr::Expr)::Union{Symbol,Expr}
        if expr.head === :(::) || expr.head === :kw
            return sign2call(expr.args[1])
        else # :parameters
            return Expr(expr.head, map(sign2call, expr.args)...)
        end # if ||
    end # function sign2call
    # find function definition
    funcdef = exprs[end]
    funcnode = findexpr(funcdef, :function)
    funcsign = findexpr(funcnode, :call)
    funcname = funcsign.args[1]
    hyfuncvar = gensym(funcname)
    callexpr::Expr = sign2call(funcsign)
    callexpr.args[1] = hyfuncvar
    # generate code
    return esc(
        quote
            let $(exprs[1:end-1]...)
                $funcdef
                global const $hyfuncvar::typeof($funcname) = $funcname
            end # let $vars
            $(funcnode.args[1]) = $callexpr
        end # return quote
    ) # esc
end # macro persistent

# Progress operations
function display_time(time::Float64)::String
    if isfinite(time) # remaining time unknown
        timeint = round(Int, time)
        min = fld(timeint, 60)
        sec = timeint % 60
        return string(min, ':', string(sec; pad=2))
    else # !isfinite(time)
        return "-:--"
    end # if isfinite, else
end # function display_time

function output!(prog::Progress, feedargs::Tuple=())::Nothing
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
    println(SS.styled"{bold,region,warning:$(prog.title)}")
    prog.lines += 1 # !
    # get bar and info strings
    elapsed = display_time(now-prog.started)
    if prog.current >= prog.total # done
        isdone = true
        # progress
        barstr = SS.annotatedstring(
            # current/total
            lpad(SS.styled"{success:$(prog.current)}", ndigits(prog.total) + 1), '/', prog.total,
            # bar
            " [", SS.styled"""{bold,success:$(repeat("━", prog.barwidth))}""", "] ",
            # percentage
            lpad(SS.styled"{success:$(round(Int, prog.current/prog.total*100))%}", 5)
        ) # SS.annotatedstring
        # time and speed
        speed = prog.current / (now-prog.started)
        togo = display_time((prog.total-prog.current) / speed)
        prompt = SS.styled"{success:{bold:Done} ✓}"
    else # in progress
        # progress
        done = floor(Int, prog.current / prog.total * prog.barwidth) # number of chars to fill =
        barstr = SS.annotatedstring(
            # current/total
            lpad(SS.styled"{info:$(prog.current)}", ndigits(prog.total) + 1), '/', prog.total,
            # bar
            " [",
            SS.styled"""{info:{bold:$(repeat("━", done))}❯}""",
            SS.styled"""{note:$(repeat("─", max(prog.barwidth-done-1, 0)))}""",
            "] ",
            # percentage
            lpad(SS.styled"{info:$(round(prog.current/prog.total*100; digits=1))%}", 5)
        ) # SS.annotatedstring
        # time and speed
        speed = (prog.current-prog.last) / (now-prog.updated)
        togo = display_time((prog.total-prog.current) / speed)
        prompt = SS.styled"{info:{bold:In progress} $(prog.runners[mod1(prog.updates, 4)])}"
    end # if >=, else
    prog.last = prog.current # !
    prog.updated = now # !
    prog.updates += 1 # !
    if !isfinite(speed) # no speed info
        spdstr = "-/sec"
    elseif (speed >= 1.0) || (iszero(speed)) # speed > 1.0
        spdstr = string(round(speed; digits=2), "/sec")
    else # speed < 1.0
        spdstr = string(round(1.0/speed; digits=2), "sec/1")
    end # if >, elseif, else
    timespeed = SS.annotatedstring(
        ' ',
        SS.styled"{$(isdone ? :success : :info):$elapsed}", "/", SS.styled"{note:-$togo}",
        ' ',
        spdstr
    ) # SS.annotatedstring
    infopaddings = repeat(" ", max(prog.width-length(timespeed)-length(prompt), 1))
    # output bar and info
    println(barstr)
    prog.lines += 1 # !
    println(timespeed, infopaddings, prompt)
    prog.lines += 1 # !
    # update user custom info
    userstr::String = prog.infofeed(feedargs...)
    userstrvec = split(userstr, '\n')
    annotatedvec = map((s -> SS.styled" {note:$s}"), userstrvec)
    foreach(s -> println(s), annotatedvec)
    prog.lines += length(annotatedvec) # !
    return nothing
end # function output

function update!(prog::Progress, current::Int=prog.current+1; feedargs::Tuple=())::Nothing
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

# Safehouse operations
function safehouse(modu::Module=Main, name::Symbol=:SAFEHOUSE)::Safehouse{modu}
    if isdefined(modu, name)
        existed = getproperty(modu, name)
        if existed isa Safehouse{modu} # exists and correct type
            return existed
        else # exists but not a Safehouse{modu}
            @warn "A variable named `$name` already exists in module `$modu` but is not a Safehouse. This variable has been housed in a new Safehouse with the given name `$name`."
            tempname = gensym(name) # protect existing variable
            safehouse = Safehouse{modu}(tempname)
            house!(name, safehouse) # house the existing variable
            @eval modu $name = $safehouse # overwrite existing variable
            return safehouse
        end # if isa, else
    else # create new safehouse
        return Safehouse{modu}(name)
    end # if isdefined, else
end # function safehouse

function house!(var::Symbol, safehouse::Safehouse{M}=safehouse())::Refugee{M} where M
    refugee = Refugee{M}(var)
    id = refugee.id
    (var in keys(safehouse.variables)) ? push!(safehouse.variables[var], id) : safehouse.variables[var] = [id] # !
    safehouse.refugees[id] = refugee # !
    return refugee
end # function house!

function note!(id::UInt32, note::String, safehouse::Safehouse{M}=safehouse())::Refugee{M} where M
    refugee = safehouse.refugees[id]
    refugee.note = note
    return refugee
end # function note!

function favorite!(id::UInt32, key::Symbol, safehouse::Safehouse{M}=safehouse())::Refugee{M} where M
    refugee = safehouse.refugees[id]
    safehouse.favorites[key] = refugee.id # !
    return refugee
end # function favorite!

(retrieve(id::UInt32, safehouse::Safehouse{M}=safehouse())::Refugee{M}) where M = safehouse.refugees[id]
(retrieve(var::Symbol, safehouse::Safehouse{M}=safehouse())::Vector{Refugee{M}}) where M =
    getindex.(safehouse.refugees, safehouse.variables[var])

# miscellaneous utilities
unique_id()::UInt32 = UInt32(UUIDs.uuid1().value >> 96)
(reprhex(hex::T)::String) where T<:Unsigned = repr(hex)[3:end]

iobuffer(io::IO; sizemodifier::NTuple{2,Int}=(0, 0))::IOContext = IOContext(
    IOBuffer(),
    :limit => true,
    :displaysize => displaysize(io) .+ sizemodifier,
    :compact => true,
    :color => true
)

# mean across vectors
@inline function crossmean(vecvec::Vector{Vector{T}})::Vector{T} where T<:Number
    @boundscheck if !all(length.(vecvec) .== length(vecvec[1]))
        throw(BoundsError("All vectors must be the same length."))
    end # if !
    return map((xi -> Stats.mean([vecvec[ti][xi] for ti in eachindex(vecvec)])), eachindex(vecvec[1]))
end # function crossmean

function hemispheric_mean(vec::Vector{T}, x::Vector{T})::T where T<:Number
    int = zero(T)
    for i in 1:length(x)-1
        @inbounds int += (vec[i]+vec[i+1]) * (x[i+1]-x[i]) / 2.0
    end # for i
    return int
end # function hemispheric_mean

# conditional copy in place
@inline function condset!(to::Vector{T}, from::T, cond::Function, ref::Vector{T}=to)::Vector{T} where T
    @. to[cond(ref)] = from # !
    return to
end # function condset!

@inline (condset(to::Vector{T}, from::T, cond::Function, ref::Vector{T}=to)::Vector{T}) where T =
    condset!(copy(to), from, cond, ref)

# replace entries with zeros in ref with zeros in place in v
@inline (zeroref!(v::Vector{T}, ref::Vector{T})::Vector{T}) where T<:Number = condset!(v, zero(T), iszero, ref)

end # module Utilities


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


module MIZ # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import NonlinearSolve as NlSol

# solar radiation absorbed on ice and water
@inline (
    solar!(base::Vector{T}, x::Vec, t::Float64, ::Val{true}, par::Collection{Float64})::Vector{T}
) where T<:Number = @. (base += par.ai * (par.S0 - par.S1 * x * cos(2.0*pi * t) - par.S2 * x^2))
@inline (
    solar!(base::Vector{T}, x::Vec, t::Float64, ::Val{false}, par::Collection{Float64})::Vector{T}
) where T<:Number = @. (base += (par.a0 - par.a2 * x^2) * (par.S0 - par.S1 * x * cos(2.0*pi * t) - par.S2 * x^2))

@inline solar(x::Vec, t::Float64, ice::Bool, par::Collection{Float64})::Vec = solar!(
    zeros(Float64, length(x)), x, t, Val(ice), par
)

# temperatures
@inline function Tbar!(Ti::Vector{T}, Tw::Vec, phi::Vec)::Vector{T} where T<:Number
    Ti .*= phi # !
    @. Ti += (1 - phi) * Tw # !
    return Ti
end # function Tbar!
@inline Tbar(Ti::Vec, Tw::Vec, phi::Vec)::Vec = Tbar!(copy(Ti), Tw, phi)
const T̄ = Tbar
const T̄! = Tbar!

@inline water_temp(Ew::Vec, phi::Vec, par::Collection{Float64})::Vec = @. par.Tm + Ew / ((1-phi) * par.cw)
@inline (ice_temp(T0::Vector{T}, par::Collection{Float64})::Vector{T}) where T<:Number = min.(T0, par.Tm)

function T0eq(
    T0::Vector{T},
    args::@NamedTuple{
        x::Vec, t::Float64, h::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime{F}, par::Collection{Float64}
    }
)::Vector{T} where {T<:Number, F} # T0eq
    vec = @. args.par.k * (args.par.Tm - T0) / args.h # SCM
    solar!(vec, args.x, args.t, Val(true), args.par) # solar on ice
    @. vec += (-args.par.A) - args.par.B * (T0 - args.par.Tm) # OLR
    diffusion!(vec, Tbar!(ice_temp(T0, args.par), args.Tw, args.phi), args.st, args.par) # diffusion
    vec .+= args.f # forcing
    return vec
end # function T0eq

@persistent T0::Vec=zeros(Float64, 100) function solveTi(
    x::Vec, t::Float64, h::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime{F}, par::Collection{Float64};
    verbose::Bool=false
)::Vec where F
    hp = condset(h, par.hmin, iszero) # avoid division by zero when solving T0
    if length(T0) != length(x)
        T0 = zeros(Float64, length(x)) # initialise T0
    end # if !=
    T0sol = NlSol.solve(
        NlSol.NonlinearProblem(T0eq, T0, (; x, t, h=hp, Tw, phi, f, st, par)),
        NlSol.TrustRegion();
        reltol=1e-6,
        abstol=1e-8
    )
    if !NlSol.SciMLBase.successful_retcode(T0sol) && verbose
        @warn "Solving for T0 failed at t=$t. Maximum residual $(maximum(abs.(T0sol.resid)))."
    end # if &&
    T0 = T0sol.u
    Ti = ice_temp(T0, par)
    zeroref!(Ti, h) # set Ti to 0 where no ice
    return Ti
end # function solveTi

# lateral melt rate
wlat(Tw::Vec, par::Collection{Float64})::Vec = @. par.m1 * (Tw - par.Tm^par.m2)

# concentration
function concentration(Ei::Vec, h::Vec, par::Collection{Float64})::Vec
    phi = @. -Ei / (par.Lf * h)
    zeroref!(phi, h)
    condset!(phi, 1.0, >(1.0)) # correct concentration
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
    x::Vec, t::Float64, ice::Bool, Ti::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime{F}, par::Collection{Float64}
)::Vec where F
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
    condset!(Ql, 0.0, isone, phi) # fix rounding errors
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
    t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime{F}, par::Collection{Float64};
    debug::Union{Expr,Nothing}=nothing, verbose::Bool=false
)::Collection{Vec} where F
    # update temperature
    vars.Tw = water_temp(vars.Ew, vars.phi, par) # !
    condset!(vars.Tw, 0.0, isnan) # eliminate NaNs for calculations
    vars.Ti = solveTi(st.x, t, vars.h, vars.Tw, vars.phi, f, st, par; verbose=verbose) # !
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
    if !isnothing(debug)
        vars.debug = eval(debug) # !
    end # if !=
    # set NaNs to no existence
    condset!(vars.Ti, NaN, iszero, vars.Ei)
    condset!(vars.Tw, NaN, >(0.99), vars.phi)
    return vars
end # function step

for xfunc in (identity, sin)
    precompile(solveTi, (Vec, Float64, Vec, Vec, Vec, Float64, SpaceTime{xfunc}, Collection{Float64}))
    precompile(step!, (Float64, Float64, Collection{Vec}, SpaceTime{xfunc}, Collection{Float64}))
end # for xfunc

end # module MIZ


module Classic # EnergyBalanceModel.

using ..Utilities, ..Infrastructure

import LinearAlgebra as LA, SparseArrays as SA

@persistent(
    cg_tau::Float64, dt_tau::Float64, dc::Float64, kappa::Matrix{Float64},
    S::Matrix{Float64}, M::Float64, aw::Vec, kLf::Float64,
    id::UInt64 = UInt64(0),

    @inline function get_statics(st::SpaceTime{F}, par::Collection{Float64})::@NamedTuple{
        cg_tau::Float64, dt_tau::Float64, dc::Float64, kappa::Matrix{Float64},
        S::Matrix{Float64}, M::Float64, aw::Vec, kLf::Float64
    } where F
        if id != hash((st, par)) # recompute only if st or par changed
            # Difinitions for implicit scheme for Tg
            cg_tau = par.cg / par.tau
            dt_tau = st.dt / par.tau
            dc = dt_tau * cg_tau
            kappa = (1+dt_tau) * LA.I(st.nx) - st.dt * par.D * get_diffop(st.nx) / par.cg
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
) # @persistent

function step!(
    t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime{F}, par::Collection{Float64};
    debug::Union{Expr,Nothing}=nothing
)::Collection{Vec} where F
    # get static variables
    stat = get_statics(st, par)
    # get time index
    i = round(Int, mod1((t + st.dt/2.0) * st.nt, st.nt))
    # forcing
    alpha = @. stat.aw * (vars.E>0.0) + par.ai * (vars.E<0.0) # WE15 Eq. (4)
    C = @. alpha*stat.S[:,i] + stat.cg_tau*vars.Tg - par.A + f
    # surface temperature
    T0 = @. C / (stat.M - stat.kLf/vars.E) # WE15 Eq. (A3)
    vars.T = @. vars.E/par.cw * (vars.E>=0) + T0 * (vars.E<0.0)*(T0<0.0) # WE15 Eq. (9)
    # Forward Euler for E
    @. vars.E += st.dt * (C - stat.M*vars.T + par.Fb) # WE15 Eq. (A2)
    # Implicit Euler for Tg
    vars.Tg =
        (stat.kappa - SA.spdiagm(stat.dc ./ (stat.M .- stat.kLf./vars.E) .* (T0.<0.0).*(vars.E.<0.0))) \
        (
            vars.Tg +
            (
                stat.dt_tau * (vars.E/par.cw.*(vars.E.>=0) +
                (par.ai*view(stat.S, :, i+1) .- par.A .+ f) ./ (stat.M .- stat.kLf./vars.E) .* (T0.<0.0).*(vars.E.<0.0))
            )
        ) # () # vars.Tg # WE15 Eq. (A1)
    # Infer ice thickness
    vars.h = @. -vars.E / par.Lf * (vars.E<0.0)
    # debug
    if !isnothing(debug)
        vars.debug = eval(debug) # !
    end # if !=
    return vars
end # function step

precompile(step!, (Float64, Float64, Collection{Vec}, SpaceTime{identity}, Collection{Float64}))

end # module ClassicEBM


module IO # EnergyBalanceModel.

using ..Utilities, ..Infrastructure

import Dates, JLD2, Makie, TimeZones as TZ

export save, load!

function unsafesave(sols, path::String; spwarn::Bool=false)::String
    if !spwarn
        @warn "`unsafesave` may overwrite existing files. Use `save` instead."
    end # if !
    JLD2.save_object(path, sols)
    return path
end

function unsafesave(plt::Makie.Figure, path::String; spwarn::Bool=false, kwargs...)::String
    if !spwarn
        @warn "`unsafesave` may overwrite existing files. Use `save` instead."
    end # if !
    Makie.save(path, plt; kwargs...)
    return path
end

function save(obj, path::String=joinpath(pwd(), string(reprhex(unique_id()), ".dat")))::String
    if isfile(path)
        modified = Dates.format(
            TZ.astimezone(
                TZ.ZonedDateTime(Dates.unix2datetime(mtime(path)), TZ.tz"UTC"),
                TZ.localzone()
            ),
            Dates.dateformat"on d u Y at HH:MM:SS"
        ) # Dates.format
        nameext = splitext(path)
        newpath = string(nameext[1], '_', reprhex(unique_id()), nameext[2])
        @warn "File $path already exists. Last modified $modified. The EXISTING file has been renamed to $newpath."
        mv(path, newpath)
    end # if isfile
    return unsafesave(obj, path; spwarn=true)
end # function save

function unsafeload(path::String; spwarn::Bool=false)
    if !spwarn
        @warn "`unsafeload` could overwrite existing variables. Use `load!` instead."
    end # if !
    return JLD2.load_object(path)
end # function unsafeload

function load!(to::Symbol, path::String, modu::Module=Main; house::Symbol=:SAFEHOUSE)
    if isdefined(modu, to)
        refugee = house!(to, safehouse(modu, house))
        @warn "Variable `$to` already defined in $modu. The existing value has been stored in safehouse `$modu.$safehouse` with ID $(reprhex(refugee.id))."
    end # if isdefined
    loaded = unsafeload(path; spwarn=true)
    @eval modu $to = $loaded
    return loaded
end # function load!

end # module IO


module Plot # EnergyBalanceModel.

using ..Utilities, ..Infrastructure

import Makie

export Layout, backend
export plot_raw, plot_avg, plot_seasonal

struct Layout{T}
    vars::Matrix{T}
    titles::Matrix{AbstractString}
    function Layout{T}(vars::Matrix{T}, titles::Matrix{AbstractString}) where T
        if size(vars) != size(titles)
            throw(ArgumentError("Size of vars and titles must be the same."))
        end # if !=
        return new{T}(vars, titles)
    end # function Layout
    Layout(vars::Matrix{T}, titles::Matrix{AbstractString}) where T = Layout{T}(vars, titles)
end

(Base.size(layout::Layout{T})::NTuple{2,Int}) where T = size(layout.vars)
(Base.axes(layout::Layout{T}, dim::Int)::Base.OneTo{Int64}) where T = axes(layout.vars, dim)
(Base.eachindex(layout::Layout{T})::Base.OneTo{Int64}) where T = eachindex(layout.vars)
(Base.getindex(layout::Layout{T}, inx...)::@NamedTuple{var::T, title::AbstractString}) where T =
    (var=layout.vars[inx...], title=layout.titles[inx...])

const miz_layout = Layout(
    [
        :Ew :Ei :E
        :Tw :Ti :T
        :h  :D  :phi
    ],
    AbstractString[
        Makie.L"$E_w$ ($\mathrm{J\,m^{-2}}$)"  Makie.L"$E_i$ ($\mathrm{J\,m^{-2}}$)"       Makie.L"$E$ ($\mathrm{J\,m^{-2}}$)"
        Makie.L"$T_w$ ($\mathrm{\degree\!C}$)" Makie.L"$T_i$ ($\mathrm{\degree\!C}$)"      Makie.L"$T$ ($\mathrm{\degree\!C}$)"
        Makie.L"$\bar{h}$ ($\mathrm{m}$)"      Makie.L"$\bar{\mathcal{D}}$ ($\mathrm{m}$)" Makie.L"\varphi"
    ]
)

const classic_layout = Layout(
    [:E :T :h],
    AbstractString[Makie.L"$E$ ($\mathrm{J\,m^{-2}}$)" Makie.L"$T$ ($\mathrm{\degree\!C}$)" Makie.L"$h$ ($\mathrm{m}$)"]
)

function init_backend(backend::Symbol=:GLMakie)::Module
    if !isdefined(Plot, backend)
        @eval import $backend
    end # if !
    modu = eval(backend)
    if Makie.current_backend() !== modu
        backend===:GLMakie ? modu.activate!(; focus_on_show=true) : modu.activate!()
    end # if !==
    return modu
end # function init_backend

backend()::Union{Module,Missing} = Makie.current_backend()
backend(backend::Symbol)::Module = init_backend(backend)

function contourf_tiles(t::Vector{T}, x::Vec, layout::Layout{Matrix{Float64}})::Makie.Figure where T<:Real
    fig = Makie.Figure()
    for row in axes(layout, 1), col in axes(layout, 2)
        subfig = fig[row,col]
        ax = Makie.Axis(
            subfig[1,1];
            title=layout[row,col].title,
            xlabel=(row==lastindex(layout, 1) ? Makie.L"$t$ ($\mathrm{y}$)" : ""),
            ylabel=(col==1 ? Makie.L"x" : ""),
            limits=(nothing, (0, 1))
        )
        ctr = Makie.contourf!(ax, t, x, layout[row,col].var)
        Makie.Colorbar(subfig[1,2], ctr)
    end # for row, col
    @eval Main innerfig = $fig
    return fig
end # function contourf_tiles

matricify(vecvec::Vector{Vec})::Matrix{Float64} = permutedims(reduce(hcat, vecvec))

function plot_raw(
    sols::Solutions{F,C};
    backend::Symbol=:GLMakie,
    layout::Layout{Symbol}=(:phi in keys(sols.raw) ? miz_layout : classic_layout)
)::Makie.Figure where {F, C}
    init_backend(backend)
    datatitle = Layout(Matrix{Matrix{Float64}}(undef, size(layout)), layout.titles)
    @simd for inx in eachindex(layout)
        datatitle.vars[inx] = matricify(getproperty(sols.raw, layout[inx].var))
    end # for inx
    return contourf_tiles(sols.ts, sols.spacetime.x, datatitle)
end # function plot_raw

function plot_avg(
    sols::Solutions{F,C};
    backend::Symbol=:GLMakie,
    layout::Layout{Symbol}=(:phi in keys(sols.raw) ? miz_layout : classic_layout)
)::Makie.Figure where {F, C}
    init_backend(backend)
    datatitle = Layout(Matrix{Matrix{Float64}}(undef, size(layout)), layout.titles)
    @simd for inx in eachindex(layout)
        datatitle.vars[inx] = matricify(getproperty(sols.seasonal.avg, layout[inx].var))
    end # for inx
    return contourf_tiles(collect(1:sols.spacetime.dur), sols.spacetime.x, datatitle)
end # function plot_avg

function plot_seasonal(
    sols::Solutions{F,false};
    xfunc::Function=(
        (sols::Solutions{F,false}, year::Int) ->
            hemispheric_mean(sols.seasonal.avg.T[year], sols.spacetime.x)
    ),
    yfunc::Function=(
            :phi in keys(sols.raw) ?
                (
                    (sols::Solutions{F,false}, season::Symbol, year::Int) ->
                        2.0*pi * hemispheric_mean(getproperty(sols.seasonal, season).phi[year], sols.spacetime.x)
                ) :
                (
                    (sols, season, year) ->
                        2.0*pi * hemispheric_mean((getproperty(sols.seasonal, season).E[year]<0.0), sols.spacetime.x)
                ) # ? :
        ), # () # yfunc
    title::AbstractString="Ice covered area",
    xlabel::AbstractString=Makie.L"$\tilde{\mathsf{T}}$ ($\mathrm{\degree\!C}$)",
    ylabel::AbstractString=Makie.L"A_i$",
    backend::Symbol=:GLMakie
)::Makie.Figure where F
    init_backend(backend)
    xdata = xfunc.(Ref(sols), sols.spacetime.dur)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1,1]; title=title, xlabel=xlabel, ylabel=ylabel)
    groups = (
        Warming=Vector{Makie.Lines{Tuple{Vector{Makie.Point{2,Float64}}}}}(),
        Cooling=Vector{Makie.Lines{Tuple{Vector{Makie.Point{2,Float64}}}}}()
    )
    for (domain, group, inx, colour) in zip(
        keys(groups),
        values(groups),
        (sols.forcing.domain[2]:sols.forcing.domain[3], sols.forcing.domain[4]:sols.forcing.domain[5]),
        (Makie.wong_colors()[6], Makie.wong_colors()[1])
    ), season in (:avg, :winter, :summer)
        width = 1.0
        if season === :avg
            width += domain===:Warming ? 2.0 : 1.0
        end # if ===
        push!(
            group,
            Makie.lines!(
                ax, xdata[inx], yfunc.(Ref(sols), Ref(season), sols.spacetime.dur[inx]);
                color=colour, linewidth=width, linestyle=(season===:summer ? :dash : :solid)
            )
        ) # push!
    end # for domain, inx, colour, season
    Makie.Legend(
        fig[1,2], values(groups), fill(["mean", "winter", "summer"], 2), keys(groups)
    )
    return fig
end # function plot_seasonal

end # module Plot


using .Utilities, .Infrastructure, .IO, .Plot

export Vec, Collection, SpaceTime, Forcing, Solutions
export miz_paramset, classic_paramset, default_parameters
export integrate
export safehouse, save, load!
export Layout, backend, plot_raw, plot_avg, plot_seasonal

end # module EnergyBalanceModel
