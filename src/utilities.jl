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
