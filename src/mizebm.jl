
module MIZEBM # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import LinearAlgebra as LA
import SparseArrays as SA

# solar radiation absorbed on ice and water
solar(x::Vec, t::Float64, ::Val{:ice}, par::Collection{Float64})::Vec =
    @. par.ai * (par.S0 - par.S1 * x * cos(2.0pi * t) - par.S2 * x^2.0)
solar(x::Vec, t::Float64, ::Val{:water}, par::Collection{Float64})::Vec =
    @. (par.a0 - par.a2 * x^2.0) * (par.S0 - par.S1 * x * cos(2.0pi * t) - par.S2 * x^2.0)
solar(x::Vec, t::Float64, surface::Symbol, par::Collection{Float64})::Vec = solar(x, t, Val(surface), par)

# phi-weighted average
weighted_avg(vi::Vec, vw::Vec, phi::Vec)::Vec = @. vi*phi + (1-phi)vw

# temperatures
water_temp(Ew::Vec, phi::Vec, par::Collection{Float64})::Vec = @. par.Tm + Ew / ((1-phi)par.cw)
water_temp_nonan(Ew::Vec, phi::Vec, par::Collection{Float64})::Vec = condset!(
    water_temp(Ew, phi, par), 0.0, ==(1.0), phi
)

ice_temp(T0::Vec, par::Collection{Float64})::Vec = min.(T0, par.Tm)

solveT0(x::Vec, t::Float64, h::Vec, Tg::Vec, Tw::Vec, phi::Vec, f::Float64, par::Collection{Float64})::Vec =
    @. (par.Tm * (par.B + par.k/h) + par.cg/par.tau * (Tg - (1-phi)Tw) + $(solar(x, t, :ice, par)) - par.A + f) /
    (par.B + par.k/h + par.cg/par.tau * phi)

function stepTg!(
    t::Float64, Tg::Vec, h::Vec, T0::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime{F}, par::Collection{Float64}
)::Vec where F
    frez = @. (T0<par.Tm) & (h>0.0)
    watr = .~frez
    M = par.B*LA.I + par.k*LA.I * SA.spdiagm(inv.(h)) + par.cg/par.tau * SA.spdiagm(phi)
    Tg .= (
            (1+st.dt/par.tau)LA.I -
            st.dt*par.D/par.cg * get_diffop(st) -
            (st.dt*par.cg/par.tau^2.0 * SA.spdiagm(phi) / M)SA.spdiagm(frez)
        ) \ (
            Tg +
            st.dt/par.tau * (
                (LA.I-SA.spdiagm(phi))Tw +
                (
                    SA.spdiagm(phi) * (
                        par.Tm * (par.B*LA.I + par.k*LA.I * SA.spdiagm(inv.(h))) -
                        par.cg/par.tau * (LA.I-SA.spdiagm(phi))SA.spdiagm(Tw) +
                        SA.spdiagm(solar(st.x, t, :ice, par)) - par.A*LA.I + f*LA.I
                    ) / M
                )frez +
                par.Tm * SA.spdiagm(phi)*watr
            )
        )
    return Tg
end # function stepTg!

# lateral melt rate
wlat(Tw::Vec, par::Collection{Float64})::Vec = @. par.m1 * (Tw - par.Tm^par.m2)

# concentration
function concentration(Ei::Vec, h::Vec, par::Collection{Float64})::Vec
    phi = @. -Ei / (par.Lf * h)
    zeroref!(phi, h)
    condset!(phi, 1.0, >(1.0)) # correct concentration
    # phi = @. Float64(Ei<0.0) # reproducing WE15
end # function concentration

# floe number
function num(D::Vec, phi::Vec, par::Collection{Float64})::Vec
    n = @. phi / (par.alpha * D^2.0)
    zeroref!(n, D)
    return n
end # function num

# lead region area
function area_lead(D::Vec, phi::Vec, n::Vec, par::Collection{Float64})::Vec
    ring = @. par.alpha * n * ((D + 2.0par.rl)^2.0 - D^2.0)
    return min.(ring, 1.0 .- phi)
end # function area_lead

# fluxes
function vert_flux(
    t::Float64, surface::Symbol, Tg::Vec, Tbar::Vec, f::Float64, st::SpaceTime{F}, par::Collection{Float64}
)::Vec where F
    L = @. par.A + par.B * (Tbar - par.Tm) # OLR
    return solar(st.x, t, surface, par) .- L .+ par.cg/par.tau * (Tg-Tbar) .+ par.Fb .+ f
end # function vert_flux

function lat_flux(h::Vec, D::Vec, Tw::Vec, phi::Vec, par::Collection{Float64})::Vec
    Flat = @. phi * h * par.Lf * $(wlat(Tw, par)) * pi / (par.alpha*D)
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
    Ql = @. Al / (1.0-phi) * psiEw
    condset!(Ql, 0.0, isone, phi) # fix rounding errors
    Qp = psiEw - Ql
    return (; Ql, Qp)
end # function split_psiEw

psinplus(Qp::Vec, par::Collection{Float64})::Vec = -Qp / (par.Lf * par.alpha * par.Dmin^2.0 * par.hmin)

function average(f::Vec, fn::Float64, n::Vec, dn::Vec)::Vec
    total = n .+ dn
    avgd = @. (n*f + dn*fn) / total
    zeroref!(avgd, total)
    # return f # reproducing WE15
    return avgd
end # function average

# differential equations
Ei_t(phi::Vec, Fvi::Vec, Flat::Vec)::Vec = @. phi * Fvi + Flat
Ew_t(phi::Vec, Fvw::Vec, Flat::Vec)::Vec = @. (1.0-phi)Fvw - Flat
h_t(Fvi::Vec, par::Collection{Float64})::Vec = -1.0/par.Lf * Fvi
function D_t(h::Vec, D::Vec, Tw::Vec, phi::Vec, Ql::Vec, par::Collection{Float64})::Vec
    lat_melt = -pi / 2.0 * par.alpha * wlat(Tw, par)
    lat_grow = @. -D / (2.0 * par.Lf * h * phi) * Ql
    weld = @. par.kappa * par.alpha / 4.0 * phi * D^3.0
    zeroref!(lat_grow, h)
    return @. lat_melt + lat_grow + weld
end # function D_t

function Infrastructure.initialise(
    ::MIZModel, st::SpaceTime{F}, forcing::Forcing{C}, par::Collection{Float64}, init::Collection{Vec};
    lastonly::Bool=true
)::Tuple{Collection{Vec}, Solutions{MIZModel,F,C}, Solutions{MIZModel,F,C}} where {F, C}
    # create storages
    vars = deepcopy(init)
    solvars = Set{Symbol}((:Ei, :Ew, :D, :h, :E, :Ti, :Tw, :T, :phi, :n))
    sols = Solutions{MIZModel}(st, forcing, par, init, solvars, lastonly) # final output
    annusol = Solutions{MIZModel}(st, forcing, par, init, solvars, true) # for annual means (internal use)
    # compute phi and Tw
    vars.nextphi = concentration(vars.Ei, vars.h, par)
    vars.nextTw = water_temp(vars.Ew, vars.nextphi, par)
    vars.nextT0 = solveT0(st.x, st.T[1], vars.h, vars.Tg, vars.nextTw, vars.nextphi, forcing(st.T[1]), par)
    condset!(vars.nextTw, 0.0, isnan) # eliminate NaNs for calculations
    return (vars, sols, annusol)
end # function initialise

forward_euler(var::Vec, grad::Vec, dt::Float64)::Vec = @. var + grad*dt

function Infrastructure.step!(
    ::MIZModel, t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime{F}, par::Collection{Float64}
)::Collection{Vec} where F
    # copy next variables to current
    vars.phi = copy(vars.nextphi)
    vars.Tw = copy(vars.nextTw)
    T0 = copy(vars.nextT0)
    # compute diagnostic variables
    vars.Ti = ice_temp(T0, par)
    condset!(vars.Ti, 0.0, isnan) # eliminate NaNs for calculations
    vars.T = weighted_avg(vars.Ti, vars.Tw, vars.phi)
    vars.n = num(vars.D, vars.phi, par)
    # calculate fluxes
    Fvi = vert_flux(t, :ice, vars.Tg, vars.T, f, st, par)
    Fvw = vert_flux(t, :water, vars.Tg, vars.T, f, st, par)
    Flat = lat_flux(vars.h, vars.D, vars.Tw, vars.phi, par)
    # update enthalpy
    rEi = forward_euler(vars.Ei, Ei_t(vars.phi, Fvi, Flat), st.dt)
    rEw = forward_euler(vars.Ew, Ew_t(vars.phi, Fvw, Flat), st.dt)
    Epsidt = redistributeE(rEi, rEw)
    vars.Ei = Epsidt.Ei # !
    vars.Ew = Epsidt.Ew # !
    vars.E = weighted_avg(vars.Ei, vars.Ew, vars.phi) # !
    # update floe size and thickness
    Al = area_lead(vars.D, vars.phi, vars.n, par)
    Qlp = split_psiEw(Epsidt.psiEwdt/st.dt, vars.phi, Al)
    dn = st.dt * psinplus(Qlp.Qp, par) # number of new pancakes
    lasth = copy(vars.h) # save for D
    vars.h = forward_euler(
        average(vars.h, par.hmin, vars.n, dn), # new pancakes
        h_t(Fvi, par),
        st.dt
    ) # !
    vars.D = forward_euler(
        average(vars.D, par.Dmin, vars.n, dn), # new pancakes
        D_t(lasth, vars.D, vars.Tw, vars.phi, Qlp.Ql, par),
        st.dt
    ) # !
    clamp!(vars.h, 0.0, Inf) # avoid overshooting to negative thickness
    zeroref!(vars.h, vars.Ei) # restrict non-existence
    clamp!(vars.D, par.Dmin, par.Dmax)
    zeroref!(vars.D, vars.Ei) # restrict non-existence
    # update variables for Tg
    vars.nextphi = concentration(vars.Ei, vars.h, par) # !
    vars.nextTw = water_temp(vars.Ew, vars.nextphi, par) # !
    condset!(vars.nextTw, 0.0, !isfinite) # eliminate NaNs for calculations
    vars.nextT0 = solveT0(st.x, t, vars.h, vars.Tg, vars.nextTw, vars.nextphi, f, par)
    vars.Tg = stepTg!(t, vars.Tg, vars.h, vars.nextT0, vars.nextTw, vars.nextphi, f, st, par) # !
    # set NaNs to no existence
    condset!(vars.Ti, NaN, iszero, vars.Ei)
    condset!(vars.Tw, NaN, >(0.99), vars.phi)
    return vars
end # function Infrastructure.step!

precompile(solveT0, (Vec, Float64, Vec, Vec, Vec, Vec, Float64, Collection{Float64}))
for xfunc in (identity, sin)
    precompile(stepTg!, (Float64, Vec, Vec, Vec, Vec, Vec, Float64, SpaceTime{xfunc}, Collection{Float64}))
    precompile(
        Infrastructure.step!,
        (MIZModel, Float64, Float64, Collection{Vec}, SpaceTime{xfunc}, Collection{Float64})
    )
end # for xfunc

end # module MIZEBM
