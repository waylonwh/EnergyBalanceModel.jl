
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
