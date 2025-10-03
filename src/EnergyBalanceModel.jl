module EnergyBalanceModel # .


module Utilities # EnergyBalanceModel.

import StyledStrings, Statistics

export Progress, update!
export inxmean, condcopy!, zeronan!, zeroref!

# progress bar
mutable struct Progress
    title::String
    total::Int
    current::Int
    last::Int # last printed progress
    started::Float64 # start time
    updated::Float64 # last external update time
    freq::Float64 # external update frequency
    infofeed::Function # Function(done::Bool, args...)::AbstractString
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
        ) # new(
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
    end # if isfinite(time), else
    return str
end # function display_time

function update!(prog::Progress, current::Int=prog.current + 1, feedargs::Tuple{Vararg{Any}}=())::Nothing
    # internal update
    prog.current = current
    # initialise if not started
    if isnan(prog.started)
        prog.started = time()
        prog.updated = time() - prog.freq # force immediate external update
    end # if isnan(prog.started)
    # external update
    isdone = false
    if (prog.current >= prog.total) || (time() - prog.updated >= prog.freq)
        now = time()
        # clear previous lines
        while prog.lines > 0
            print("\033[A\033[2K") # move up one line and clear the line
            prog.lines -= 1
        end # while prog.lines > 0
        if prog.current > prog.total
            return nothing # avoid over-updating
        end # if prog.current > prog.total
        # title
        println(StyledStrings.styled"{bold,region,warning:$(prog.title)}")
        prog.lines += 1
        # get bar and info strings
        elapsed = display_time(now - prog.started)
        if prog.current == prog.total # done
            isdone = true
            # progress
            barstr = StyledStrings.annotatedstring(
                # current/total
                lpad(StyledStrings.styled"{success:$(prog.current)}", ndigits(prog.total) + 1), '/', prog.total,
                # bar
                " [", StyledStrings.styled"""{success:$(repeat("=", prog.barwidth))}""", "] ",
                # percentage
                lpad(StyledStrings.styled"{success:$(round(Int, prog.current/prog.total*100))%}", 5)
            ) # StyledStrings.annotatedstring(
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
            ) # StyledStrings.annotatedstring(
            # time and speed
            speed = (prog.current - prog.last) / (now - prog.updated)
            togo = display_time((prog.total - prog.current) / speed)
            prompt = StyledStrings.styled"{info:{bold:In progress} $(prog.runners[mod1(prog.updates, 4)])}"
        end # if prog.current == prog.total, else
        prog.last = prog.current
        prog.updated = now
        prog.updates += 1
        if !isfinite(speed) # no speed info
            spdstr = "-/sec"
        elseif (speed >= 1.0) || (iszero(speed)) # speed > 1.0
            spdstr = string(round(speed, digits=2), "/sec")
        else # speed < 1.0
            spdstr = string(round(1.0 / speed, digits=2), "sec/1")
        end # if speed > 1.0, elseif, else
        timespeed = StyledStrings.annotatedstring(
            ' ',
            StyledStrings.styled"{$(isdone ? :success : :info):$elapsed}",
            "/",
            StyledStrings.styled"{note:-$togo}",
            ' ',
            spdstr
        ) # StyledStrings.annotatedstring(
        infopaddings = repeat(" ", max(prog.width - length(timespeed) - length(prompt), 1))
        # output bar and info
        println(barstr)
        prog.lines += 1
        println(timespeed, infopaddings, prompt)
        prog.lines += 1
    end # if (prog.current == prog.total) || (now-prog.updated >= prog.freq)
    # update user custom info
    userstr = prog.infofeed(isdone, feedargs...)
    userstrvec = split(userstr)
    annotatedvec = map((s -> StyledStrings.styled" {note:$s}"), userstrvec)
    foreach(s::AbstractString -> println(s), annotatedvec)
    prog.lines += length(annotatedvec)
    return nothing
end # function update!

@inline function inxmean(vecvec::Vector{Vector{T}})::Vector{T} where T<:Number
    @boundscheck if !all(length.(vecvec) .== length(vecvec[1]))
        throw(BoundsError("All vectors must be the same length."))
    end # if !all(length.(vecvec) .== length(vecvec[1]))
    return  map((xi::Int -> Statistics.mean([vecvec[ti][xi] for ti in eachindex(vecvec)])), eachindex(vecvec[1]))
end # function inxmean

# conditional copy in place
@inline function condcopy!(to::Vector{T}, from::T, cond::Function, ref::Vector{T}=to)::Vector{T} where T
    @. to[cond(ref)] = from
    return to
end # function condcopy!

@inline function condcopy!(to::Vector{T}, from::Vector{T}, cond::Function, ref::Vector{T}=to)::Vector{T} where T
    @. to[cond(ref)] = from[cond(ref)]
    return to
end # function condcopy!

# replace NaNs with zeros in place
@inline (zeronan!(v::Vector{T}, ref::Vector{T}=v)::Vector{T}) where T<:Number = condcopy!(v, zero(T), isnan, ref)
# replace entries with zeros in ref with zeros in place in v
@inline (zeroref!(v::Vector{T}, ref::Vector{T})::Vector{T}) where T<:Number = condcopy!(v, zero(T), iszero, ref)

end # module Utilities


module Infrastructure # EnergyBalanceModel.

using ..Utilities

import SparseArrays

export Floats, SymbSet, Vec, Vecs, SpaceTime, Solutions, Forcing
export default_params, miz_paramset, classic_paramset
export get_defaultpar, diffusion, D∇², annual_mean

const SymbSet = Set{Symbol}
const Vec = Vector{Float64}

macro asstruct(def::Expr)
    # parse definition
    name = def.args[1] # =
    dicttype = def.args[2] # =
    valtype = def.args[2].args[3] # =.curly
    # generate code
    return esc(
        quote
            struct $name
                dict::$dicttype
                $name(args...) = new($dicttype(args...))
            end # struct $name

            Base.getproperty(obj::$name, key::Symbol)::$valtype = getindex(getfield(obj, :dict), key)
            Base.setproperty!(obj::$name, key::Symbol, val::$valtype)::$dicttype = setindex!(
                getfield(obj, :dict), val, key
            )
            Base.keys(obj::$name)::AbstractSet{Symbol} = keys(getfield(obj, :dict))
        end # quote
    ) # esc(
end # macro asstruct

@asstruct Floats = Dict{Symbol,Float64}
@asstruct Vecs = Dict{Symbol,Vec}
@asstruct NestedVecs = Dict{Symbol,Vector{Vec}} # solution storage unit

struct SpaceTime
    nx::Int # number of evenly spaced latitudinal gridboxes (equator to pole)
    dx::Float64 # grid box width
    x::Vec # grid
    dur::Int # duration of simulation in years
    nt::Int # number of timesteps per year (limited by numerical stability)
    dt::Float64 # timestep
    t::Vec # time vector in a year
    T::AbstractRange{Float64} # time vector for entire simulation
    winter::@NamedTuple{t::Float64,inx::Int}
    summer::@NamedTuple{t::Float64,inx::Int}
    diffop::SparseArrays.SparseMatrixCSC{Float64,Int64} # difference operator matrix

    function get_diffop(nx::Int, dx::Float64)::SparseArrays.SparseMatrixCSC{Float64,Int64}
        xb = dx : dx: 1.0-dx
        lambda = @. (1 - xb^2) / dx^2
        l1 = pushfirst!(-copy(lambda), 0.0)
        l2 = push!(-copy(lambda), 0.0)
        l3 = -l1 - l2
        return SparseArrays.spdiagm(-1 => -l1[2:nx], 0 => -l3, 1 => -l2[1:nx-1])
    end # function get_diffop

    function SpaceTime(nx::Int, nt::Int, dur::Int, winter::Float64=0.26125, summer::Float64=0.77375)
        dx = 1.0 / nx
        x = collect(range(dx/2.0, 1.0 - dx/2.0, nx))
        dt = 1.0 / nt
        t = collect(range(dt/2.0, 1.0 - dt/2.0, nt))
        T = dt/2.0 : dt : dur - dt/2.0
        winterinx = round(Int, nt*winter)
        summerinx = round(Int, nt*summer)
        diffop = get_diffop(nx, dx)
        return new(nx, dx, x, dur, nt, dt, t, T, (t=winter, inx=winterinx), (t=summer, inx=summerinx), diffop)
    end # function SpaceTime
end # struct SpaceTime

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
        end # if T < forcing.domain[2], elseif*3, else
    end # if forcing.flat, else
    return f
end # function (forcing::Forcing)

struct Solutions
    spacetime::SpaceTime # space and time which solutions are defined on
    ts::Vec # time vector for stored solution
    forcing::Forcing # climate forcing
    parameters::Floats # model parameters
    initconds::Vecs # initial conditions
    lastonly::Bool # store only last year of solution
    debug::Expr # store debug variables
    raw::NestedVecs # solution storage
    seasonal::@NamedTuple{winter::NestedVecs, summer::NestedVecs, avg::NestedVecs} # seasonal peak and annual avg

    function Solutions(
        st::SpaceTime, forcing::Forcing, par::Floats, init::Vecs, vars::SymbSet,
        lastonly::Bool=true;
        debug::Expr=Expr(:block)
    ) # Solutions(
        if lastonly
            dur_store = 1
            ts = st.dur-1.0 + st.dt/2.0 : st.dt : st.dur - st.dt/2.0
        else # !lastonly
            dur_store = st.dur
            ts = st.dt/2.0 : st.dt : st.dur - st.dt/2.0
        end # if lastonly, else
        if debug != Expr(:block)
            push!(vars, :debug)
        end # if debug != Expr(:block)
        # construct raw solution storage
        solraw = NestedVecs()
        foreach((var::Symbol -> setproperty!(solraw, var, Vector{Vec}(undef, length(ts)))), vars)
        # construct seasonal solution storage template
        seasonaltemp = NestedVecs()
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
        ) # new(
    end # function Solutions
end # struct Solutions

# default parameter values
const default_params = Floats(
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
    :kappa => 0.01  * 31536000 # floe welding parameter
) # Floats(

# Floats used in each model
const miz_paramset = SymbSet(
    [
        :D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :Tm, :m1, :m2, :alpha,
        :rl, :Dmin, :Dmax, :hmin, :kappa
    ]
)
const classic_paramset = SymbSet(
    [:D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :F, :cg, :tau]
)

# Create a parameter dictionary from default values for a given Set
function get_defaultpar(paramset::SymbSet)::Floats
    setvec = collect(paramset)
    return Floats(setvec .=> getproperty.(Ref(default_params), setvec))
end # function get_defaultpar

# calculate diffusion operator matrix
@inline (diffusion(f::VT, st::SpaceTime, par::Floats)::VT) where VT<:Vector{<:Number} = par.D * st.diffop * f
const D∇² = diffusion

# calculate annual mean
function annual_mean(annusol::Solutions)::Vecs
    # calculate annual mean for each variable except temperatures
    means = Vecs()
    foreach(
        (var::Symbol -> setproperty!(means, var, inxmean(getproperty(annusol.raw, var)))),
        keys(annusol.raw)
    )
    return means
end # function annual_mean

annual_mean(forcing::Forcing, st::SpaceTime, year::Int)::Float64 = Statistics.mean(forcing.(year-1 .+ st.t))

end # module Infrastructure


module MIZEBM # EnergyBalanceModel.

module Components # EnergyBalanceModel.MIZEBM.

using ...Infrastructure, ...Utilities

import NonlinearSolve

export insolation, coalbedo
export water_temp, solveTi, Tbar, T̄
export wlat, vert_flux, lat_flux, redistributeE
export split_psiEw, psinplus, area_lead, average
export concentration, num
export Ei_t, Ew_t, h_t, D_t

# Functions of space and time
@inline insolation(x::Vec, t::Float64, par::Floats)::Vec = @. par.S0 - par.S1 * x * cos(2 * pi * t) - par.S2 * x^2
@inline coalbedo(x::Vec, ice::Bool, par::Floats)::Vec = ice ? fill(par.ai, length(x)) : @. par.a0 - par.a2 * x^2

# temperatures
@inline function Tbar(Ti::VT, Tw::Vec, phi::Vec)::VT where VT<:Vector{<:Number}
    zeronan!(Ti)
    zeronan!(Tw)
    return @. phi * Ti + (1 - phi) * Tw
end # function Tbar
const T̄ = Tbar

water_temp(Ew::Vec, phi::Vec, par::Floats)::Vec = @. par.Tm + Ew / ((1-phi) * par.cw)
@inline (ice_temp(T0::VT, par::Floats)::VT) where VT<:Vector{<:Number} = min.(T0, par.Tm)

function T0eq(
    T0::VT,
    args::@NamedTuple{
        x::Vec, t::Float64, h::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime, par::Floats
    }
)::VT where VT<:Vector{<:Number} # T0eq(
    conduction = @. args.par.k * (args.par.Tm - T0) / args.h
    solar = coalbedo(args.x, true, args.par) .* insolation(args.x, args.t, args.par)
    olr = @. -args.par.A - args.par.B * (T0 - args.par.Tm)
    dlap = diffusion(Tbar(ice_temp(T0, args.par), args.Tw, args.phi), args.st, args.par)
    forcing = args.f
    return @. conduction + solar + olr + dlap + forcing
end # function T0eq

solveTi =
    let T0::Vec = zeros(Float64, 100) # let T0 be a persistent variable
        function solveTi(x::Vec, t::Float64, h::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime, par::Floats)::Vec
            iced = (h .> 0.0)
            h_pos = copy(h)
            @. h_pos[!iced] = 1.0 # avoid division by zero
            zeronan!(Tw) # avoid NaNs in Tbar
            if length(T0) != length(x)
                T0 = zeros(length(x)) # initialise T0
            end # if length(T0) != length(x)
            T0sol = NonlinearSolve.solve(
                NonlinearSolve.NonlinearProblem(T0eq, T0, (; x, t, h=h_pos, Tw, phi, f, st, par)),
                NonlinearSolve.TrustRegion()
            )
            # TODO test for solving failure
            T0 = T0sol.u
            Ti = ice_temp(T0, par)
            @. Ti[!iced] = NaN # reset to NaN where no ice
            return Ti
        end # function solveTi
    end # let T0

# lateral melt rate
function wlat(Tw::Vec, par::Floats)::Vec
    w = @. par.m1 * (Tw - par.Tm^par.m2)
    zeronan!(w, Tw)
    return w
end # function wlat

# concentration
function concentration(Ei::Vec, h::Vec, par::Floats)::Vec
    phi = @. -Ei / (par.Lf * h)
    zeroref!(phi, h)
    condcopy!(phi, 1.0, >(1.0)) # correct concentration
    # phi = map((e -> Float64(e<0.0)), Ei) # reproducing WE15
    return phi
end # function concentration

# floe number
function num(D::Vec, phi::Vec, par::Floats)::Vec
    n = @. phi / (par.alpha * D^2)
    zeroref!(n, D)
    return n
end # function num

# lead region area
function area_lead(D::Vec, phi::Vec, n::Vec, par::Floats)::Vec
    ring = @. par.alpha * n * ((D + 2*par.rl)^2 - D^2)
    return min.(ring, 1.0.-phi)
end # function area_lead

# fluxes
function vert_flux(
    x::Vec, t::Float64, ice::Bool, Ti::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime, par::Floats
)::Vec # vert_flux(
    zeronan!(Ti)
    zeronan!(Tw)
    L = @. par.A + par.B * ($Tbar(Ti, Tw, phi) - par.Tm) # OLR
    return @. $coalbedo(x, ice, par)*$insolation(x, t, par) - L + $(diffusion(Tbar(Ti, Tw, phi), st, par)) + par.Fb + f
end # function vert_flux

function lat_flux(h::Vec, D::Vec, Tw::Vec, phi::Vec, par::Floats)::Vec
    Flat = @. phi * h * par.Lf * $wlat(Tw, par) * pi / (par.alpha * D)
    zeroref!(Flat, D)
    zeronan!(Flat, Tw)
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

psinplus(Qp::Vec, par::Floats)::Vec = -Qp / (par.Lf * par.alpha * par.Dmin^2 * par.hmin)

function average(f::Vec, fn::Float64, n::Vec, dn::Vec)::Vec
    total = n .+ dn
    avgd = @. (n*f + dn*fn) / total
    zeroref!(avgd, total)
    return avgd
end # function average

# differential equations
Ei_t(phi::Vec, Fvi::Vec, Flat::Vec)::Vec = @. phi * Fvi + Flat
Ew_t(phi::Vec, Fvw::Vec, Flat::Vec)::Vec = @. (1 - phi) * Fvw - Flat
h_t(Fvi::Vec, par::Floats)::Vec = -1/par.Lf * Fvi
function D_t(h::Vec, D::Vec, Tw::Vec, phi::Vec, Ql::Vec, par::Floats)::Vec
    lat_melt = -pi / (2 * par.alpha) * wlat(Tw, par)
    lat_grow = @. -D / (2 * par.Lf * h * phi) * Ql
    weld = @. par.kappa * par.alpha / 4 * phi * D^3
    zeroref!(lat_grow, h)
    return @. lat_melt + lat_grow + weld
end # function D_t

end # module Components


module Integration # EnergyBalanceModel.MIZEBM.

using ...Infrastructure, ...Utilities, ..Components

export integrate

forward_euler(var::Vec, grad::Vec, dt::Float64)::Vec = @. var + grad * dt

function step!(
    t::Float64, f::Float64, vars::Vecs, st::SpaceTime, par::Floats;
    debug::Expr=Expr(:block)
)::Vecs
    # update temperature
    vars.Tw = water_temp(vars.Ew, vars.phi, par) # !
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
    clamp!(vars.D, par.Dmin, par.Dmax) # TODO test if needed
    zeroref!(vars.D, vars.Ei) # correct round off errors
    rh = forward_euler(vars.h, h_t(Fvi, par), st.dt)
    clamp!(rh, 0.0, Inf) # avoid overshooting to negative thickness
    vars.h = average(rh, par.hmin, vars.n, dn) # new pancakes # !
    # update concentration
    vars.phi = concentration(vars.Ei, vars.h, par) # !
    vars.E = @. vars.phi * vars.Ei + (1 - vars.phi) * vars.Ew # !
    vars.T = Tbar(vars.Ti, vars.Tw, vars.phi) # !
    # debug
    if debug != Expr(:block)
        vars.debug = eval(debug) # !
    end # if debug != Expr(:block)
    return vars
end # function step

function savesol!(sols::Solutions, annusol::Solutions, vars::Vecs, tinx::Int)::Solutions
    varscp = deepcopy(vars) # avoid reference issues
    year = ceil(Int, sols.spacetime.T[tinx])
    ti = mod1(tinx, sols.spacetime.nt) # index of time in the year
    # save raw data to annual
    foreach(keys(annusol.raw)) do var::Symbol
        getproperty(annusol.raw, var)[ti] = getproperty(varscp, var)
    end # foreach() do
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
    end # if !sols.lastonly, elseif
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
    end # if ti % sols.spacetime.nt == sols.spacetime.winter.inx, elseif*2
    return sols
end # function savesol!

function integrate(
    st::SpaceTime, forcing::Forcing, par::Floats, init::Vecs;
    lastonly::Bool=true, debug::Expr=Expr(:block)
)::Solutions
    # initialise
    vars = deepcopy(init)
    solvars = SymbSet([:Ei, :Ew, :E, :Ti, :Tw, :T, :h, :D, :phi, :n])
    sols = Solutions(st, forcing, par, init, solvars, lastonly, debug=debug)
    annusol = Solutions(st, forcing, par, init, solvars, true, debug=debug) # for calculating annual means
    progress = Progress(length(st.T), "Integrating")
    update!(progress)
    # loop over time
    for ti in eachindex(st.T)
        step!(st.t[mod1(ti, st.nt)], forcing(st.T[ti]), vars, st, par, debug=debug)
        savesol!(sols, annusol, vars, ti)
        update!(progress)
    end # for ti in eachindex(st.T)
    return sols
end # function integrate

end # module Integration

import .Integration: integrate

export integrate

end # module MIZEBM


module ClassicEBM # EnergyBalanceModel.

end # module ClassicEBM


end # module EnergyBalanceModel
