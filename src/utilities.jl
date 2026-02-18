module Utilities # EnergyBalanceModel.

import Statistics as Stats, StyledStrings as SS

export Progress, update!
export @isdebugging, @persistent, iobuffer
export condset, condset!, crossmean, zeroref!

# add new function introduced in Julia 1.12
@static if VERSION < v"1.12"
    export ctruncate
    ctruncate(x, _...) = x
end # if <

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

# progress bar
mutable struct Progress
    title::String # title of the progress
    from::Float64 # starting point
    current::Float64 # current progress
    to::Float64 # end point
    isdone::Bool # whether the progress is complete
    digits::Int # number of digits to display in numbers
    var::String # variable name for nums
    started::Float64 # start time
    history::Vector{NTuple{3,Float64}} # [(time, percentage, num),]
    infofeed::Function # Function(args...)::String
    feedargs::Tuple # arguments for infofeed
    width::Int # number of characters wide, including progress texts
    barwidth::Int # width of the progress bar
    numswidth::Int # width of the number display (current/total)
    lines::Int # lines printed
    updates::Int # number of external updates
    timer::Timer # timer for automatic updates

    function Progress(
        from::Float64, to::Float64, var::String="", title::String="Progress";
        infofeed::Function=Returns(""), digits::Union{Int,Nothing}=nothing, width::Int=60
    )
        if isnothing(digits) # determine minor digits
            if all(<=(1)∘abs, (from, to)) # 0.xxx
                digits = 3
            elseif any(>(1000)∘abs, (from, to)) # 1000.x
                digits = 1
            else # 1.xx 10.xx 100.xx
                digits = 2
            end # if all, elseif any, else
        end # if isnothing
        maxnumlen = maximum(length∘string∘Base.Fix1(trunc, Int), (from, to)) + 1 + digits
        numswidth = 2maxnumlen + 1
        barwidth = width - numswidth - 4 - 5 # (var = )1.23/10.00 [=> ] xx.x%
        isempty(var) || (barwidth -= length(var) + 3)
        return new(
            #            ↓current  ↓isdone        started↓           ↓history
            title, from, from, to, false, digits, var, NaN, Vector{NTuple{3,Float64}}(),
            infofeed, (), width, barwidth, 0, 0, Timer(Returns(nothing), 0)
            #         ^feedargs       lines^  ^updates     ^tiemr
        ) # new
    end # function Progress

    Progress(total::Float64, arg...; kwargs...) = Progress((0.0, total), arg...; kwargs...)
end # struct Progress

function speed!(prog::Progress)::Float64
    length(prog.history) < 2 && return NaN # not enough history to determine speed
    # keep history within 1% of current progress
    while prog.history[2][2] < prog.history[end][2] - 0.01
        deleteat!(prog.history, 1)
    end # while <
    return (prog.history[end][3] - prog.history[1][3]) / (prog.history[end][1] - prog.history[1][1])
end # function speed!

title_anstr(prog::Progress)::Base.AnnotatedString{String} =
    SS.styled"{bold,region,warning:$(prog.title)}"

function nums_anstr(prog::Progress, isdone::Bool)::Base.AnnotatedString{String}
    varstr = isempty(prog.var) ? "" : string(prog.var, " = ")
    color = isdone ? :success : :info
    return SS.annotatedstring(
        varstr, SS.styled"{$color:$(prog.current)}", '/',
        round(prog.to; digits=prog.digits), # current stage
    ) # SS.annotatedstring
end # function nums_anstr

bar_anstr(::Val{true}, barwidth::Int, _)::Base.AnnotatedString{String} = SS.annotatedstring(
    " [", SS.styled"{bold,success:$(repeat('━', barwidth))}", "] ",
)
bar_anstr(::Val{false}, barwidth::Int, filling::Float64)::Base.AnnotatedString{String} = SS.annotatedstring(
    " [", SS.styled"{info:{bold:$(repeat('━', filling))}❯}",
    SS.styled"{note:$(repeat('─', max(barwidth-filling-1, 0)))}", "] ",
)

percentage_anstr(::Val{true}, _)::Base.AnnotatedString{String} = SS.styled"{success:100%}"
percentage_anstr(::Val{false}, percentage::Float64)::Base.AnnotatedString{String} =
    SS.styled"{info:$(round(percentage*100; digits=1))%}"

function time_str(time::Float64)::String
    if isfinite(time) # remaining time unknown
        timeint = round(Int, time)
        min = fld(timeint, 60)
        sec = timeint % 60
        return string(min, ':', string(sec; pad=2))
    else # !isfinite(time)
        return "-:--"
    end # if isfinite, else
end # function time_str

function speed_str(speed::Float64)::String
    if !isfinite(speed) # no speed info
        return "-/sec"
    elseif speed>=1 || iszero(speed) # speed > 1
        return string(round(speed; digits=2), "/sec")
    else # speed < 1
        return string(round(1/speed; digits=2), "sec/1")
    end # if >, elseif, else
end # function speed_str

function timespeed_anstr(
    ::Val{true}, prog::Progress, now::Float64, _
)::Base.AnnotatedString{String}
    elapsed = now - prog.started
    speed = (prog.to - prog.from) / elapsed
    return SS.annotatedstring(
        SS.styled"{success:$(time_str(elapsed))}", '/', SS.styled"{note:-$(time_str(0.0))}",
        ' ', speed_str(speed)
    )
end # function timespeed_anstr

function timespeed_anstr(
    ::Val{false}, prog::Progress, now::Float64, speed::Float64
)::Base.AnnotatedString{String}
    elapsed = now - prog.started
    togo = (prog.to - prog.current) / speed
    return SS.annotatedstring(
        SS.styled"{info:$(time_str(elapsed))}", '/', SS.styled"{note:-$(time_str(togo))}",
        ' ', speed_str(speed)
    )
end # function timespeed_anstr

prompt_anstr(::Val{true}, _)::Base.AnnotatedString{String} = SS.styled"{success:{bold:Done} ✓}"
function prompt_anstr(::Val{false}, updates::Int)::Base.AnnotatedString{String}
    runner = ('◓', '◑', '◒', '◐')[mod1(updates, 4)]
    return SS.styled"{info:{bold:In progress} $runner}"
end

function output!(prog::Progress)::Nothing
    now = time()
    prog.updates += 1 # !
    percentage = min((prog.current - prog.from) / (prog.to - prog.from), 1.0)
    filling = floor(Int, percentage * prog.barwidth)
    push!(prog.history, (now, percentage, prog.current)) # !
    speed = speed!(prog)
    # clear previous lines
    print(repeat("\033[A\033[2K", prog.lines)) # move up and clear lines
    prog.lines = 0 # !
    # title
    println(title_anstr(prog)) # [title]  Stage 1/1
    prog.lines += 1 # !
    # prgress & bar
    print(lpad(nums_anstr(prog, prog.isdone), prog.numswidth)) # 1.00/2.00
    print(bar_anstr(Val(prog.isdone), prog.barwidth, filling)) # [━❯─]
    println(lpad(percentage_anstr(Val(prog.isdone), percentage), 5)) # 12.3% 100%
    prog.lines += 1 # !
    # time, speed and prompt
    timespeed = timespeed_anstr(Val(prog.isdone), prog, now, speed) # 0:17/-0:26 1.23/sec
    prompt = prompt_anstr(Val(prog.isdone), prog.updates) # Done ✓  In progress ◓
    paddings = repeat(' ', max(prog.width-length(timespeed)-length(prompt), 1))
    println(timespeed, paddings, prompt)
    prog.lines += 1 # !
    # user custom info
    userstr::String = prog.infofeed(prog.feedargs...)
    userstrvec = split(userstr, '\n')
    annotatedvec = map((s -> SS.styled"{note:$s}"), userstrvec)
    foreach(println, annotatedvec)
    prog.lines += length(annotatedvec) # !
    return nothing
end # function output

function start!(prog::Progress; feedargs::Tuple=())::Nothing
    prog.started = time() # !
    prog.feedargs = feedargs # !
    prog.timer = Timer(_ -> output!(prog), 0; interval=0.2)
    return nothing
end # function start!

function update!(
    prog::Progress, current::Float64; feedargs::Tuple=()
)::Nothing
    prog.current = current # !
    prog.feedargs = feedargs # !
    prog.isdone = prog.to > prog.from ? prog.current >= prog.to : prog.current <= prog.to # !
    if prog.isdone
        isopen(prog.timer) || @warn "Progress output timer is closed before completion!"
        close(prog.timer) # stop timer
        output!(prog) # final output
    end # if &&
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
@inline function crossmean(vecvec::Vector{Vector{T}})::Vector{T} where T<:Number
    @boundscheck all(length.(vecvec) .== length(vecvec[1])) ||
        throw(BoundsError("All vectors must be the same length."))
    return map((xi -> Stats.mean(vecvec[ti][xi] for ti in eachindex(vecvec))), eachindex(vecvec[1]))
end # function crossmean

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
