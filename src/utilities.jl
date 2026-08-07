module Utilities # EnergyBalanceModel.

import Statistics as Stats, StyledStrings as SS

export Progress, update!
export @isdebugging, @persistent, iobuffer
export condset!, crossmean, zeroref!

# add new function introduced in Julia 1.12
@static if VERSION < v"1.12"
    export ctruncate
    ctruncate(x, _...) = x
end # if <

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
    runners::Tuple{Vararg{Char}} # characters to use as runners
    updates::Int # number of external updates

    function Progress(
        total::Integer,
        title::String="Progress", freq::Real=1.0;
        width::Int=50, infofeed::Function=Returns("")
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

macro persistent(exprs...) # -> Expr
    # syntax tree operations
    findexpr(_, ::Symbol)::Nothing = nothing
    function findexpr(expr::Expr, head::Symbol)::Union{Expr,Nothing}
        if expr.head === head
            return expr
        elseif isempty(expr.args)
            return nothing
        else # recursively search args
            for arg in expr.args
                funcexpr = findexpr(arg, head)
                isnothing(funcexpr) || return funcexpr
            end # for arg
            return nothing
        end # if ==
    end # function findexpr
    sign2call(expr::Symbol)::Symbol = expr
    sign2call(expr::Expr)::Union{Symbol,Expr} = (expr.head===:(::) || expr.head===:kw) ?
        sign2call(expr.args[1]) : Expr(expr.head, map(sign2call, expr.args)...)
    # find function definition
    funcdef = exprs[end]
    funcnode = findexpr(funcdef, :function)
    funcsign = findexpr(funcnode, :call)
    funcname = funcsign.args[1]
    hyfuncvar = gensym(funcname)
    callexpr = sign2call(funcsign)
    callexpr.args[1] = hyfuncvar
    # generate code
    return esc(
        quote
            let $(exprs[1:end-1]...)
                $funcdef
                global const $hyfuncvar::typeof($funcname) = $funcname
            end # let $vars
            $(funcnode.args[1]) = $callexpr
        end # quote
    ) # esc
end # macro persistent

# determine is the current file or module is being debugged
macro isdebugging() # -> Bool
    dbg_env = get(ENV, "JULIA_DEBUG", "")
    if isempty(dbg_env)
        return false
    else # !isempty
        targets = split(dbg_env, ',')
        file = splitext(splitpath(string(__source__.file))[end])[1]
        mods = split(string(__module__), '.')
        return file in targets || !isempty(intersect(targets, mods))
    end # if isempty, else
end # macro isdebugging

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
    # clear previous lines
    print(repeat("\033[A\033[2K", prog.lines)) # move up and clear lines
    prog.lines = 0 # !
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
        prompt = SS.styled"{info:{bold:In progress} $(prog.runners[mod1(prog.updates, length(prog.runners))])}"
    end # if >=, else
    prog.last = prog.current # !
    prog.updated = now # !
    prog.updates += 1 # !
    if !isfinite(speed) # no speed info
        spdstr = "-/sec"
    elseif speed>=1 || iszero(speed) # speed > 1
        spdstr = string(round(speed; digits=2), "/sec")
    else # speed < 1
        spdstr = string(round(1/speed; digits=2), "sec/1")
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
    annotatedvec = map(s -> SS.styled" {note:$s}", userstrvec)
    foreach(println, annotatedvec)
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
    ((time() - prog.updated >= prog.freq) || (prog.current == prog.total)) &&
        output!(prog, feedargs)
    return nothing
end # function update!

iobuffer(io::IO; sizemodifier::NTuple{2,Int}=(0, 0))::IOContext = IOContext(
    IOBuffer(),
    :limit => true,
    :displaysize => displaysize(io) .+ sizemodifier,
    :compact => true,
    :color => true
)

# mean across vectors
@inline function crossmean(vecvec::Vector{<:Vector}) # -> Vector{Number}
    @boundscheck all(length.(vecvec) .== length(vecvec[1])) ||
        throw(BoundsError("All vectors must be the same length."))
    return map(xi -> Stats.mean(vecvec[ti][xi] for ti in eachindex(vecvec)), eachindex(vecvec[1]))
end # function crossmean

# conditional copy in place
@inline function condset!(to::Vector, from, cond::Function, ref::Vector) # -> Vector{T}
    @. to[cond(ref)] = from # !
    return to
end # function condset!

# replace entries with zeros in ref with zeros in place in v
@inline (zeroref!(v::Vector{T}, ref::Vector{T})::Vector{T}) where T = condset!(v, zero(T), iszero, ref)

end # module Utilities
