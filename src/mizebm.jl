
module MIZEBM # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import NonlinearSolve as NlSol

# solar radiation absorbed on ice and water
@inline (
    solar!(base::Vector{T}, x::Vec, t::Float64, ::Val{true}, par::Collection{Float64})::Vector{T}
) where T<:Number = @. (base += par.ai * (par.S0 - par.S1 * x * cos(2.0pi * t) - par.S2 * x^2))
@inline (
    solar!(base::Vector{T}, x::Vec, t::Float64, ::Val{false}, par::Collection{Float64})::Vector{T}
) where T<:Number = @. (base += (par.a0 - par.a2 * x^2) * (par.S0 - par.S1 * x * cos(2.0pi * t) - par.S2 * x^2))

@inline solar(x::Vec, t::Float64, ice::Bool, par::Collection{Float64})::Vec = solar!(
    zeros(Float64, length(x)), x, t, Val(ice), par
)

# phi-weighted average
@inline function weighted_avg!(vi::Vector{T}, vw::Vec, phi::Vec)::Vector{T} where T<:Number
    @. vi *= phi # !
    @. vi += (1 - phi) * vw # !
    return vi
end # function weighted_avg!
weighted_avg(vi::Vec, vw::Vec, phi::Vec)::Vec = weighted_avg!(copy(vi), vw, phi)

@inline water_temp(Ew::Vec, phi::Vec, par::Collection{Float64})::Vec = @. par.Tm + Ew / ((1-phi) * par.cw)
@inline (ice_temp(T0::Vector{T}, par::Collection{Float64})::Vector{T}) where T<:Number = min.(T0, par.Tm)

function T0eq(
    T0::Vector{T},
    args::@NamedTuple{
        t::Float64, h::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime{F}, par::Collection{Float64}
    }
)::Vector{T} where {T<:Number, F} # T0eq
    vec = @. args.par.k * (args.par.Tm - T0) / args.h # SCM
    solar!(vec, args.st.x, args.t, Val(true), args.par) # solar on ice
    @. vec += (-args.par.A) - args.par.B * (T0 - args.par.Tm) # OLR
    diffusion!(vec, weighted_avg!(ice_temp(T0, args.par), args.Tw, args.phi), args.st, args.par) # diffusion
    vec .+= args.f # forcing
    return vec
end # function T0eq

@persistent T0::Vec=zeros(Float64, 100) function solveTi(
    t::Float64, h::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime{F}, par::Collection{Float64};
    verbose::Bool=false
)::Vec where F
    hp = condset(h, par.hmin, iszero) # avoid division by zero when solving T0
    if length(T0) != length(st.x)
        T0 = zeros(Float64, length(st.x)) # initialise T0
    end # if !=
    T0sol = NlSol.solve(
        NlSol.NonlinearProblem(T0eq, T0, (; t, h=hp, Tw, phi, f, st, par)),
        NlSol.TrustRegion();
        reltol=1e-6,
        abstol=1e-8
    )
    if verbose && !NlSol.SciMLBase.successful_retcode(T0sol)
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
    ring = @. par.alpha * n * ((D + 2.0par.rl)^2 - D^2)
    return min.(ring, 1.0.-phi)
end # function area_lead

# fluxes
function vert_flux(
    t::Float64, ice::Bool, Ti::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime{F}, par::Collection{Float64}
)::Vec where F
    L = @. par.A + par.B * ($(weighted_avg(Ti, Tw, phi)) - par.Tm) # OLR
    return @. $solar(st.x, t, ice, par) - L + $(diffusion(weighted_avg(Ti, Tw, phi), st, par)) + par.Fb + f
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
    lat_melt = -pi / 2.0 * par.alpha * wlat(Tw, par)
    lat_grow = @. -D / (2.0 * par.Lf * h * phi) * Ql
    weld = @. par.kappa * par.alpha / 4 * phi * D^3
    zeroref!(lat_grow, h)
    return @. lat_melt + lat_grow + weld
end # function D_t

function Infrastructure.initialise(
    ::MIZModel, st::SpaceTime{F}, forcing::Forcing{C}, par::Collection{Float64}, init::Collection{Vec};
    lastonly::Bool=true, debug::Union{Expr,Nothing}=nothing, verbose::Bool=false
)::Tuple{Collection{Vec}, Solutions{MIZ,F,C}, Solutions{MIZ,F,C}} where {F,C}
    # create storage
    vars = deepcopy(init)
    solvars = Set{Symbol}((:Ei, :Ew, :D, :h, :E, :Ti, :Tw, :T, :phi, :n))
    sols = Solutions{MIZ}(st, forcing, par, init, solvars, lastonly; debug)
    annusol = Solutions{MIZ}(st, forcing, par, init, solvars, true; debug) # for calculating annual means
    # complete step 0
    vars.phi = concentration(vars.Ei, vars.h, par)
    vars.Tw = water_temp(vars.Ew, vars.phi, par)
    condset!(vars.Tw, 0.0, isnan) # eliminate NaNs for calculations
    vars.Ti = solveTi(0.0, vars.h, vars.Tw, vars.phi, forcing(0.0), st, par; verbose)
    vars.n = num(vars.D, vars.phi, par)
    vars.E = weighted_avg(vars.Ei, vars.Ew, vars.phi)
    vars.T = weighted_avg(vars.Ti, vars.Tw, vars.phi)
    return (vars, sols, annusol)
end # function initialise

forward_euler(var::Vec, grad::Vec, dt::Float64)::Vec = @. var + grad * dt

function Infrastructure.step!(
    ::MIZModel, t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime{F}, par::Collection{Float64};
    debug::Union{Expr,Nothing}=nothing, verbose::Bool=false
)::Collection{Vec} where F
    condset!(vars.Tw, 0.0, isnan) # eliminate NaNs for calculations
    # calculate fluxes
    Fvi = vert_flux(t, true, vars.Ti, vars.Tw, vars.phi, f, st, par)
    Fvw = vert_flux(t, false, vars.Ti, vars.Tw, vars.phi, f, st, par)
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
    # update floe number
    vars.n = num(vars.D, vars.phi, par) # !
    # update temperature
    vars.Tw = water_temp(vars.Ew, vars.phi, par) # !
    condset!(vars.Tw, 0.0, isnan) # eliminate NaNs for calculations
    vars.Ti = solveTi(t, vars.h, vars.Tw, vars.phi, f, st, par; verbose) # !
    # update total energy and temperature
    zeroref!(vars.Ei, vars.h) # correct round off errors
    vars.E = weighted_avg(vars.Ei, vars.Ew, vars.phi) # !
    vars.T = weighted_avg(vars.Ti, vars.Tw, vars.phi) # !
    # debug
    if !isnothing(debug)
        vars.debug = eval(debug) # !
    end # if !=
    # set NaNs to no existence
    condset!(vars.Ti, NaN, iszero, vars.Ei)
    condset!(vars.Tw, NaN, >(0.99), vars.phi)
    return vars
end # function Infrastructure.step!

for xfunc in (identity, sin)
    precompile(solveTi, (Vec, Float64, Vec, Vec, Vec, Float64, SpaceTime{xfunc}, Collection{Float64}))
    precompile(
        Infrastructure.step!,
        (MIZ, Float64, Float64, Collection{Vec}, SpaceTime{xfunc}, Collection{Float64})
    )
end # for xfunc

end # module MIZEBM
