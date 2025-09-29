module EnergyBalanceModel


module Infrastructure

import SparseArrays

export Parameters, ParamSet, VarVec, InitConds, SpaceTime, Solution, Forcing
export default_params, miz_paramset, classic_paramset
export get_defaultpar, get_diffop!, diffusion, D∇², condcopy!, zeronan!, zeroref!

const Parameters = Dict{Symbol,Float64}
const ParamSet = Set{Symbol}
const VarVec = Vector{Float64}
const InitConds = Dict{Symbol,VarVec}

struct SpaceTime
    nx::Int # number of evenly spaced latitudinal gridboxes (equator to pole)
    dx::Float64 # grid box width
    x::VarVec # grid
    nt::Int # number of timesteps per year (limited by numerical stability)
    dt::Float64 # timestep
    t::VarVec # time vector in a year

    function SpaceTime(nx::Int, nt::Int)
        dx = 1.0 / nx
        x = collect(range(dx / 2.0, 1.0 - dx / 2.0, nx))
        dt = 1.0 / nt
        t = collect(range(dt / 2.0, 1.0 - dt / 2.0, nt))
        return new(nx, dx, x, nt, dt, t)
    end # function SpaceTime
end # struct SpaceTime

struct Solution
    dur::Int # duration of simulation in years
    lastonly::Bool # store only last year of solution
    T::VarVec # time vector for solution storage
    sol::Dict{Symbol,Vector{VarVec}} # solution storage

    function Solution(
        dur::Int, lastonly::Bool, st::SpaceTime, vars::Set{Symbol}; debug::Bool=false
    )
        dur_store = lastonly ? 1 : dur # for calculating the number of total timesteps
        T = st.dt:st.dt:dur_store
        debug ? push!(vars, :debug) : nothing
        sol = Dict{Symbol,Vector{VarVec}}()
        foreach(var -> (sol[var] = Vector{VarVec}(undef, length(T))), vars)
        return new(dur, lastonly, T, sol)
    end # function Solution
end # struct Solution

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
            throw(ArgumentError("Warming time must be integer. Calculated: $warming."))
        # hold at peak
        @. domainvec[4:6] += holdyrs[2]
        # time to cool
        cooling = (peak - cool) / rates[2]
        isinteger(cooling) ?
            @.(domainvec[5:6] += cooling) :
            throw(ArgumentError("Cooling time must be integer. Calculated: $cooling."))
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
const miz_paramset = ParamSet(
    [
        :D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :Tm, :m1, :m2, :alpha,
        :rl, :Dmin, :Dmax, :hmin, :kappa
    ]
)
const classic_paramset = ParamSet(
    [:D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :F, :cg, :tau]
)

# Create a parameter dictionary from default values for a given Set
function get_defaultpar(paramset::ParamSet)::Parameters
    setvec = collect(paramset)
    return Parameters(setvec .=> getindex.(Ref(default_params), setvec))
end # function get_defaultpar

# calculate diffusion operator matrix
function get_diffop!(par::Parameters, st::SpaceTime)::Nothing
    xb = st.dx : st.dx : 1.0 - st.dx
    lambda = @. par[:D]/st.dx^2 * (1 - xb^2)
    l1 = pushfirst!(-copy(lambda), 0.0)
    l2 = push!(-copy(lambda), 0.0)
    l3 = -l1 - l2
    par[:diffop] = SparseArrays.spdiagm(-1 => -l1[2:st.nx], 0 => -l3, 1 => -l2[1:st.nx-1]) # !
    return nothing
end # function get_diffop!

diffusion(f::VarVec, diffop::AbstractMatrix{Float64})::VarVec = diffop * f
const D∇² = diffusion

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

end # module Infrastructure


module MIZEBM


module Components

import NonlinearSolve

using ...Infrastructure

# Functions of space and time
solar(x::VarVec, t::Float64, par::Parameters)::VarVec = @. par[:S0] - par[:S1] * x * cos(2 * pi * t) - par[:S2] * x^2
coalbedo(x::VarVec, ice::Bool, par::Parameters)::VarVec = ice ? fill(par[:ai], length(x)) : @. par[:a0] - par[:a2] * x^2

# temperatures
water_temp(Ew::VarVec, phi::VarVec, par::Parameters)::VarVec = @. par[:Tm] + Ew / ((1-phi) * par[:cw])

let T0 = fill(0.0, 1) # let T0 be a persistent variable
    ice_temp(T0::AbstractVector, par::Parameters)::VarVec = min.(T0, par[:Tm])

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
            args.par[:k] * (args.par[:Tm] - T0) / args.h + # vertical conduction in ice
            coalbedo(args.x, true, args.par) * solar(args.x, args.t, args.par) + # solar
            (-args.par[:A]) - args.par[:B] * (T0 - args.par[:Tm]) + # OLR
            diffusion(Tbar(ice_temp(T0, args.par), args.Tw, args.phi), args.par[:diffop]) + # diffusion
            args.f; # forcing
        end # @. begin

    function ice_temp(
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
        (length(T0) != length(x)) ? (T0 = fill(0.0, length(x))) : nothing # initialise T0
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
    w = @. par[:m1] * (Tw - par[:Tm]^par[:m2])
    zeronan!(w, Tw)
    return w
end # function wlat

# concentration
function concentration(Ei::VarVec, h::VarVec, par::Parameters)::VarVec
    phi = @. -Ei / (par[:Lf] * h)
    zeroref!(phi, h)
    condcopy!(phi, 1.0, >(1.0)) # correct concentration
    # phi = map(e -> Float64(e<0.0), Ei) # reproducing WE15
    return phi
end # function concentration

# floe number
function num(D::VarVec, phi::VarVec, par::Parameters)::VarVec
    n = @. phi / (par[:alpha] * D^2)
    zeroref!(n, D)
    return n
end # function num

# lead region area
function area_lead(D::VarVec, phi::VarVec, n::VarVec, par::Parameters)::VarVec
    ring = @. par[:alpha] * n * ((D + 2*par[:rl])^2 - D^2)
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
    L = @. par[:A] + par[:B] * (Tbar(Ti, Tw, phi) - par[:Tm]) # OLR
    return @. coalbedo(x, ice, par) * solar(x, t, par) - L + diffusion(Tbar(Ti, Tw, phi), par[:diffop]) + par[:Fb] + f
end # function vert_flux

function lat_flux(h::VarVec, D::VarVec, Tw::VarVec, phi::VarVec, par::Parameters)::VarVec
    Flat = @. phi * h * par[:Lf] * wlat(Tw, par) * pi / (par[:alpha] * D)
    zeroref!(Flat, D)
    zeronan!(Flat, Tw)
    return Flat
end # function lat_flux

# redistribution functions
function split_psiEw(psiEw::VarVec, phi::VarVec, Al::VarVec)::@NamedTuple{Ql::VarVec, Qp::VarVec}
    Ql = @. Al / (1-phi) * psiEw
    condcopy!(Ql, 0.0, isone, phi) # fix rounding errors
    Qp = psiEw - Ql
    return (; Ql, Qp)
end # function split_psiEw

psinplus(Qp::VarVec, par::Parameters)::VarVec = -Qp / (par[:Lf] * par[:alpha] * par[:Dmin]^2 * par[:hmin])

# differential equations
Ei_t(phi::VarVec, Fvi::VarVec, Flat::VarVec)::VarVec = @. phi * Fvi + Flat
Ew_t(phi::VarVec, Fvw::VarVec, Flat::VarVec)::VarVec = @. (1 - phi) * Fvw - Flat
h_t(Fvi::VarVec, par::Parameters)::VarVec = -1/par[:Lf] * Fvi
function D_t(h::VarVec, D::VarVec, Tw::VarVec, phi::VarVec, Ql::VarVec, par::Parameters)::VarVec
    lat_melt = -pi / (2 * par[:alpha]) * wlat(Tw, par)
    lat_grow = @. -D / (2 * par[:Lf] * h * phi) * Ql
    weld = @. par[:kappa] * par[:alpha] / 4 * phi * D^3
    zeroref!(lat_grow, h)
    return @. lat_melt + lat_grow + weld
end # function D_t

end # module Components


module Integration

import ...Infrastructure

function step()
    return nothing
end # function step

function integrate()
    return nothing
end # function integrate

end # module Integration


end # module MIZEBM


module ClassicEBM

end # module ClassicEBM


end # module EnergyBalanceModel
