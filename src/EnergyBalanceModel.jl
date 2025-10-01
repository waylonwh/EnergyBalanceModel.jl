module EnergyBalanceModel # .


module Infrastructure # EnergyBalanceModel.

import SparseArrays, Statistics

export Parameters, SymbSet, VarVec, Variables, SpaceTime, Solutions, Forcing
export default_params, miz_paramset, classic_paramset
export get_defaultpar, get_diffop!, diffusion, D∇², annual_mean

macro asstruct(def::Expr)
    # parse definition
    name = def.args[1].args[1] # const.=
    type = def.args[1].args[2].args[3] # const.=.curly
    # generate code
    return esc(
        quote
            $def
            Base.getproperty(obj::$name, key::Symbol)::$type = getindex(obj, key)
            Base.setproperty!(obj::$name, key::Symbol, val::$type)::$type = setindex!(obj, val, key)
        end # quote
    ) # esc(
end # macro asstruct

const SymbSet = Set{Symbol}
const VarVec = Vector{Float64}
@asstruct const Parameters = Dict{Symbol,Float64}
@asstruct const Variables = Dict{Symbol,VarVec}
@asstruct const Sol = Dict{Symbol,Vector{VarVec}} # solution storage unit

struct SpaceTime
    nx::Int # number of evenly spaced latitudinal gridboxes (equator to pole)
    dx::Float64 # grid box width
    x::VarVec # grid
    dur::Int # duration of simulation in years
    nt::Int # number of timesteps per year (limited by numerical stability)
    dt::Float64 # timestep
    t::VarVec # time vector in a year
    T::AbstractRange{Float64} # time vector for entire simulation
    winter::@NamedTuple{t::Float64,inx::Int}
    summer::@NamedTuple{t::Float64,inx::Int}

    function SpaceTime(nx::Int, nt::Int, dur::Int, winter::Float64=0.26125, summer::Float64=0.77375)
        dx = 1.0 / nx
        x = collect(range(dx/2.0, 1.0 - dx/2.0, nx))
        dt = 1.0 / nt
        t = collect(range(dt/2.0, 1.0 - dt/2.0, nt))
        T = dt/2.0 : dt : dur - dt/2.0
        winterinx = round(Int, nt*winter)
        summerinx = round(Int, nt*summer)
        return new(nx, dx, x, dur, nt, dt, t, T, (t=winter, inx=winterinx), (t=summer, inx=summerinx))
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
    Forcing(base::Float64) =
        new(true, base, base, base, (Inf, 0), (0.0, 0.0), (0, Inf, Inf, Inf, Inf, Inf))
    # warming/cooling forcing
    function Forcing(
        base::Float64,
        peak::Float64,
        cool::Float64,
        holdyrs::Tuple{Int,Int},
        rates::Tuple{Float64,Float64}
    ) # Forcing(
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
    ts::VarVec # time vector for stored solution
    forcing::Forcing # climate forcing
    parameters::Parameters # model parameters
    initconds::Variables # initial conditions
    lastonly::Bool # store only last year of solution
    debug::Expr # store debug variables
    raw::Sol # solution storage
    seasonal::@NamedTuple{winter::Sol, summer::Sol, avg::Sol} # average and seasonal solutions

    # create empty solution storage with given variables and size
    function emptysol(vars::SymbSet, size::Int)::Sol
        soltemp = Sol()
        foreach((var::Symbol -> setproperty!(soltemp, var, Vector{VarVec}(undef, size))), vars)
        return soltemp
    end # function emptysol

    function Solutions(
        st::SpaceTime,
        forcing::Forcing,
        par::Parameters,
        init::Variables,
        vars::SymbSet,
        lastonly::Bool=true;
        debug::Expr=:nothing
    ) # Solutions(
        if lastonly
            dur_store = 1
            ts = st.dur-1.0 + st.dt/2.0 : st.dt : st.dur - st.dt/2.0
        else # !lastonly
            dur_store = st.dur
            ts = st.dt/2.0 : st.dt : st.dur - st.dt/2.0
        end # if lastonly, else
        if debug !== :nothing
            push!(vars, :debug)
        end
        # construct raw solution storage
        solraw = emptysol(vars, length(ts))
        # construct seasonal solution storage template
        seasonaltemp = emptysol(vars, dur_store)
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
                avg=deepcopy(seasonaltemp),
                summer=deepcopy(seasonaltemp),
                winter=deepcopy(seasonaltemp)
            ) # ( # seasonal
        ) # new(
    end # function Solutions
end # struct Solutions

# default parameter values
const default_params = Parameters(
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
) # Parameters(

# Parameters used in each model
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
function get_defaultpar(paramset::SymbSet)::Parameters
    setvec = collect(paramset)
    return Parameters(setvec .=> getindex.(Ref(default_params), setvec))
end # function get_defaultpar

# calculate diffusion operator matrix
function get_diffop!(par::Parameters, st::SpaceTime)::Nothing
    xb = st.dx : st.dx : 1.0 - st.dx
    lambda = @. par.D/st.dx^2 * (1 - xb^2)
    l1 = pushfirst!(-copy(lambda), 0.0)
    l2 = push!(-copy(lambda), 0.0)
    l3 = -l1 - l2
    par.diffop = SparseArrays.spdiagm(-1 => -l1[2:st.nx], 0 => -l3, 1 => -l2[1:st.nx-1]) # !
    return nothing
end # function get_diffop!

diffusion(f::VarVec, diffop::AbstractMatrix{Float64})::VarVec = diffop * f
const D∇² = diffusion

# calculate annual mean
function annual_mean(annusol::Solutions)::Variables
    # calculate annual mean for each variable except temperatures
    means = Variables()
    foreach((var::Symbol -> setproperty!(means, var, Statistics.mean.(zip(annusol.raw.var...)))), keys(annusol.raw))
    return means
end

annual_mean(forcing::Forcing, st::SpaceTime, year::Int)::Float64 = Statistics.mean(forcing.(year-1 .+ st.t))


end # module Infrastructure


module Utilities # EnergyBalanceModel.

import StyledStrings

export Progress, update!
export condcopy!, zeronan!, zeroref!

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
        title::String="Progress",
        freq::Float64=1.0,
        infofeed::Function=(done::Bool->"");
        width::Int=25
    )
        barwidth = width - (ndigits(total)*2 + 1) - 2 - 5 - 3 # current/total [=> ] xx.x%
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
    timeint = round(Int, time)
    min = fld(timeint, 60)
    sec = timeint % 60
    return string(min, ':', string(sec, pad=2))
end # function display_time

function update!(prog::Progress, current::Int=prog.current+1, feedargs::Tuple{Vararg{Any}}=())::Nothing
    # internal update
    prog.current = current
    # initialise if not started
    if prog.started === NaN
        prog.started = time()
        prog.updated = time() - prog.freq # force immediate external update
    end
    # external update
    if (prog.current >= prog.total) || (time()-prog.updated >= prog.freq)
        now = time()
        # clear previous lines
        while prog.lines > 0
            print("\033[A\033[2K") # move up one line and clear the line
            prog.lines -= 1
        end # while prog.lines > 0
        if prog.current > prog.total
            return nothing # avoid over-updating
        end
        # title
        println(StyledStrings.styled"{bold,region,warning:$prog.title}")
        prog.lines += 1
        # get bar and info strings
        elapsed = display_time(now - prog.started)
        prog.last = prog.current
        prog.updated = now
        prog.updates += 1
        if prog.current == prog.total # done
            isdone = true
            # progress
            barstr = StyledStrings.annotatedstring(
                # current/total
                lpad(StyledStrings.styled"{success:$prog.current}", ndigits(prog.total)+1),
                '/',
                prog.total,
                # bar
                " [",
                StyledStrings.styled"{success:$('='^prog.barwidth)}",
                "] ",
                # percentage
                lpad(
                    StyledStrings.styled"{success:$(round(Int, prog.current/prog.total*100))%}",
                    5
                )
            ) # StyledStrings.annotatedstring(
            # time and speed
            speed = prog.current / elapsed
            togo = display_time((prog.total - prog.current) / speed)
            prompt = StyledStrings.styled"{success:{bold:Done} ✓}"
        else # in progress
            isdone = false
            # progress
            done = floor(Int, prog.current/prog.total * prog.barwidth) # number of chars to fill =
            barstr = StyledStrings.annotatedstring(
                # current/total
                lpad(StyledStrings.styled"{info:$prog.current}", ndigits(prog.total)+1),
                '/',
                prog.total,
                # bar
                " [",
                StyledStrings.styled"{info:$('='^done)>}",
                " "^(prog.barwidth-done-1),
                "] ",
                # percentage
                lpad(
                    StyledStrings.styled"{info:$(round(prog.current/prog.total*100, digits=1))%}",
                    5
                )
            ) # StyledStrings.annotatedstring(
            # time and speed
            speed = (prog.current-prog.last) / (now-prog.updated)
            togo = display_time((prog.total-prog.current) / speed)
            prompt = StyledStrings.styled"{info:{bold:In progress} $prog.runners[mod1(prog.updates, 4)]}"
        end # if prog.current == prog.total, else
        timespeed = StyledStrings.annotatedstring(
            ' ',
            StyledStrings.styled"{$(isdone ? :success : :info):$elapsed}",
            "/",
            StyledStrings.styled"{note:-$togo}",
            ' ',
            speed > 1.0 ?
                string(round(speed, digits=2), "/sec") :
                string(round(1.0 / speed, digits=2), "sec/1")
        ) # StyledStrings.annotatedstring(
        infopaddings = " "^(prog.width - length(timespeed) - length(prompt))
        # output bar and info
        println(barstr)
        prog.lines += 1
        println(timespeed, infopaddings, prompt)
        prog.lines += 1
    end # if (prog.current == prog.total) || (now-prog.updated >= prog.freq)
    # update user custom info
    userstr = prog.infofeed(isdone, feedargs...)
    userstrvec = split(userstr)
    annotatedvec = map((s::AbstractString -> StyledStrings.styled" {note:$s}"), userstrvec)
    foreach(s::AbstractString -> println(s), annotatedvec)
    prog.lines += length(annotatedvec)
    return nothing
end # function update!

# conditional copy in place
function condcopy!(to::Vector{T}, from::T, cond::Function, ref::Vector{T}=to)::Vector{T} where {T<:Number}
    @. to[cond(ref)] = from
    return to
end # function condcopy!

function condcopy!(to::Vector{T}, from::Vector{T}, cond::Function, ref::Vector{T}=to)::Vector{T} where {T<:Number}
    @. to[cond(ref)] = from[cond(ref)]
    return to
end # function condcopy!

# replace NaNs with zeros in place
(zeronan!(v::Vector{T}, ref::Vector{T}=v)::Vector{T}) where {T<:Number} = condcopy!(v, 0.0, isnan, ref)
# replace entries with zeros in ref with zeros in place in v
(zeroref!(v::Vector{T}, ref::Vector{T})::Vector{T}) where {T<:Number} = condcopy!(v, 0.0, iszero, ref)

end # module Utilities


module MIZEBM # EnergyBalanceModel.


module Components # EnergyBalanceModel.MIZEBM.

using ...Infrastructure, ...Utilities

import NonlinearSolve

export solar, coalbedo
export water_temp, ice_temp, Tbar, T̄
export wlat, vert_flux, lat_flux, redistributeE
export split_psiEw, psinplus, area_lead, average
export concentration, num
export Ei_t, Ew_t, h_t, D_t

# Functions of space and time
solar(x::VarVec, t::Float64, par::Parameters)::VarVec = @. par.S0 - par.S1 * x * cos(2 * pi * t) - par.S2 * x^2
coalbedo(x::VarVec, ice::Bool, par::Parameters)::VarVec = ice ? fill(par.ai, length(x)) : @. par.a0 - par.a2 * x^2

# temperatures
water_temp(Ew::VarVec, phi::VarVec, par::Parameters)::VarVec = @. par.Tm + Ew / ((1-phi) * par.cw)

let T0 = fill(0.0, 1) # let T0 be a persistent variable
    ice_temp(T0::AbstractVector, par::Parameters)::VarVec = min.(T0, par.Tm)

    T0eq(
        T0::AbstractVector,
        args::@NamedTuple{
            x::VarVec,
            t::Float64,
            h::VarVec,
            Tw::VarVec,
            phi::VarVec,
            f::Float64,
            par::Parameters
        } # @NamedTuple{
    )::AbstractVector =
        @. begin # T0eq(
            args.par.k * (args.par.Tm - T0) / args.h + # vertical conduction in ice
            coalbedo(args.x, true, args.par) * solar(args.x, args.t, args.par) + # solar
            (-args.par.A) - args.par.B * (T0 - args.par.Tm) + # OLR
            diffusion(Tbar(ice_temp(T0, args.par), args.Tw, args.phi), args.par.diffop) + # diffusion
            args.f; # forcing
        end # @. begin

    global function ice_temp(
        x::VarVec,
        t::Float64,
        h::VarVec,
        Tw::VarVec,
        phi::VarVec,
        f::Float64,
        par::Parameters
    )::VarVec # ice_temp(
        iced = (h .> 0.0)
        h_pos = copy(h)
        @. h_pos[!iced] = 1.0 # avoid division by zero
        zeronan!(Tw) # avoid NaNs in Tbar
        if length(T0) != length(x)
            T0 = fill(0.0, length(x))
        end
        T0sol = NonlinearSolve.solve(
            NonlinearSolve.NonlinearProblem(T0eq, T0, (; x, t, h=h_pos, Tw, phi, f, par))
        )
        # TODO test for solving failure
        T0 = T0sol.u
        Ti = ice_temp(T0, par)
        @. Ti[!iced] = NaN # reset to NaN where no ice
        return Ti
    end # function ice_temp
end # let T0

function Tbar(Ti::VarVec, Tw::VarVec, phi::VarVec)::VarVec
    zeronan!(Ti)
    zeronan!(Tw)
    return @. phi * Ti + (1 - phi) * Tw
end # function Tbar
const T̄ = Tbar

# lateral melt rate
function wlat(Tw::VarVec, par::Parameters)::VarVec
    w = @. par.m1 * (Tw - par.Tm^par.m2)
    zeronan!(w, Tw)
    return w
end # function wlat

# concentration
function concentration(Ei::VarVec, h::VarVec, par::Parameters)::VarVec
    phi = @. -Ei / (par.Lf * h)
    zeroref!(phi, h)
    condcopy!(phi, 1.0, >(1.0)) # correct concentration
    # phi = map((e -> Float64(e<0.0)), Ei) # reproducing WE15
    return phi
end # function concentration

# floe number
function num(D::VarVec, phi::VarVec, par::Parameters)::VarVec
    n = @. phi / (par.alpha * D^2)
    zeroref!(n, D)
    return n
end # function num

# lead region area
function area_lead(D::VarVec, phi::VarVec, n::VarVec, par::Parameters)::VarVec
    ring = @. par.alpha * n * ((D + 2*par.rl)^2 - D^2)
    return min.(ring, 1.0-phi)
end # function area_lead

# fluxes
function vert_flux(
    x::VarVec,
    t::Float64,
    ice::Bool,
    Ti::VarVec,
    Tw::VarVec,
    phi::VarVec,
    f::Float64,
    par::Parameters
)::VarVec # vert_flux(
    zeronan!(Ti)
    zeronan!(Tw)
    L = @. par.A + par.B * (Tbar(Ti, Tw, phi) - par.Tm) # OLR
    return @. coalbedo(x, ice, par) * solar(x, t, par) - L + diffusion(Tbar(Ti, Tw, phi), par.diffop) + par.Fb + f
end # function vert_flux

function lat_flux(h::VarVec, D::VarVec, Tw::VarVec, phi::VarVec, par::Parameters)::VarVec
    Flat = @. phi * h * par.Lf * wlat(Tw, par) * pi / (par.alpha * D)
    zeroref!(Flat, D)
    zeronan!(Flat, Tw)
    return Flat
end # function lat_flux

function redistributeE(
    rEi::VarVec, rEw::VarVec
)::@NamedTuple{Ei::VarVec, Ew::VarVec, psiEidt::VarVec, psiEwdt::VarVec}
    cEi = clamp.(rEi, -Inf, 0.0)
    cEw = clamp.(rEw, 0.0, Inf)
    psiEidt = rEi .- cEi # +
    psiEwdt = rEw .- cEw # -
    Ei = cEi .+ psiEwdt # -
    Ew = cEw .+ psiEidt # +
    return (; Ei, Ew, psiEidt, psiEwdt)
end # function redistributeE

# redistribution functions
function split_psiEw(psiEw::VarVec, phi::VarVec, Al::VarVec)::@NamedTuple{Ql::VarVec, Qp::VarVec}
    Ql = @. Al / (1-phi) * psiEw
    condcopy!(Ql, 0.0, isone, phi) # fix rounding errors
    Qp = psiEw - Ql
    return (; Ql, Qp)
end # function split_psiEw

psinplus(Qp::VarVec, par::Parameters)::VarVec = -Qp / (par.Lf * par.alpha * par.Dmin^2 * par.hmin)

function average(f::VarVec, fn::Float64, n::VarVec, dn::VarVec)::VarVec
    total = n .+ dn
    avgd = @. (n*f + dn*fn) / total
    zeroref!(avgd, total)
    return avgd
end # function average

# differential equations
Ei_t(phi::VarVec, Fvi::VarVec, Flat::VarVec)::VarVec = @. phi * Fvi + Flat
Ew_t(phi::VarVec, Fvw::VarVec, Flat::VarVec)::VarVec = @. (1 - phi) * Fvw - Flat
h_t(Fvi::VarVec, par::Parameters)::VarVec = -1/par.Lf * Fvi
function D_t(h::VarVec, D::VarVec, Tw::VarVec, phi::VarVec, Ql::VarVec, par::Parameters)::VarVec
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

forward_euler(var::VarVec, grad::VarVec, dt::Float64)::VarVec = @. var + grad * dt

function step!(
    t::Float64, f::Float64, vars::Variables, st::SpaceTime, par::Parameters;
    debug::Expr=:nothing
)::Variables
    # update temperature
    vars.Tw = water_temp(vars.Ew, vars.phi, par)
    vars.Ti = ice_temp(st.x, t, vars.h, vars.Tw, vars.phi, f, par)
    # update floe number
    vars.n = num(vars.D, vars.phi, par)
    # calculate fluxes
    Fvi = vert_flux(st.x, t, true, vars.Ti, vars.Tw, vars.phi, f, par)
    Fvw = vert_flux(st.x, t, false, vars.Ti, vars.Tw, vars.phi, f, par)
    Flat = lat_flux(vars.h, vars.D, vars.Tw, vars.phi, par)
    # update enthalpy
    rEi = forward_euler(vars.Ei, Ei_t(vars.phi, Fvi, Flat), st.dt)
    rEw = forward_euler(vars.Ew, Ew_t(vars.phi, Fvw, Flat), st.dt)
    Epsidt = redistributeE(rEi, rEw)
    vars.Ei = Epsidt.Ei
    vars.Ew = Epsidt.Ew
    # update floe size and thickness
    Al = area_lead(vars.D, vars.phi, vars.n, par)
    Qlp = split_psiEw(Epsidt.psiEwdt/st.dt, vars.phi, Al)
    dn = st.dt * psinplus(Qlp.Qp, par) # number of new pancakes
    rD = forward_euler(vars.D, D_t(vars.h, vars.D, vars.Tw, vars.phi, Qlp.Ql, par), st.dt)
    vars.D = average(rD, par.Dmin, vars.n, dn) # new pancakes
    # clamp!.(vars.D, par.Dmin, par.Dmax) # TODO test if needed
    # zeroref!(vars.D, vars.Ei) # TODO test if needed
    rh = forward_euler(vars.h, h_t(Fvi, par), st.dt)
    clamp!.(rh, 0.0, Inf) # avoid overshooting to negative thickness
    vars.h = average(rh, par.hmin, vars.n, dn) # new pancakes
    # update concentration
    vars.phi = concentration(vars.Ei, vars.h, par)
    vars.E = @. vars.phi * vars.Ei + (1 - vars.phi) * vars.Ew
    vars.T = Tbar(vars.Ti, vars.Tw, vars.phi)
    # debug
    if debug !== :nothing
        vars.debug = eval(debug)
    end
    return vars
end # function step

function savesol!(sols::Solutions, annusol::Solutions, vars::Variables, ti::Int)::Solutions
    varscopy = deepcopy(vars) # avoid reference issues
    # save raw data to annual
    foreach((var::Symbol -> setindex!(annusol.raw, varscopy.var, mod1(ti, sols.spacetime.nt))), keys(annusol.raw))
    # save raw data
    if !sols.lastonly # save all raw data
        foreach((var::Symbol -> setindex!(sols.raw, varscopy.var, ti)), keys(sols.raw))
    elseif ti > length(sols.spacetime.T) - sols.spacetime.nt # save the raw data of the last year
        inx = ti - (length(sols.spacetime.T) - sols.spacetime.nt)
        foreach((var::Symbol -> setindex!(sols.raw, varscopy.var, inx)), keys(sols.raw))
    end # if !sols.lastonly, elseif
    # save seasonal data
    if ti % sols.spacetime.nt == sols.spacetime.winter.inx
        foreach(
            (var::Symbol -> setindex!(sols.seasonal.winter, varscopy.var, ceil(Int, ti))),
            keys(sols.seasonal.winter)
        )
    elseif ti % sols.spacetime.nt == sols.spacetime.summer.inx
        foreach(
            (var::Symbol -> setindex!(sols.seasonal.summer, varscopy.var, ceil(Int, ti))),
            keys(sols.seasonal.summer)
        )
    elseif ti % sols.spacetime.nt == 0 # calculate annual average
        means = annual_mean(annusol)
        foreach(
            (
                var::Symbol ->
                    setindex!(sols.seasonal.avg, means.var, ceil(Int, sols.spacetime.T[ti]))
            ),
            keys(sols.seasonal.avg)
        ) # foreach(
    end # if ti % sols.spacetime.nt == sols.spacetime.winter.inx, elseif*2
    return sols
end # function savesol!

function integrate(
    st::SpaceTime, forcing::Forcing, par::Parameters, init::Variables;
    lastonly::Bool=true, debug::Expr=:nothing
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


end # module MIZEBM


module ClassicEBM # EnergyBalanceModel.

end # module ClassicEBM


end # module EnergyBalanceModel
