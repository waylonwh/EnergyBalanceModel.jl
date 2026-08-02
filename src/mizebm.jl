module MIZEBM # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import LinearAlgebra as LA
import SparseArrays as SA

# solar radiation absorbed on ice and water
solar(x::Vec, t::Float64, ::Val{:ice}, par::Par)::Vec =
    @. par.ai * (par.S0 - par.S1 * x * cos(2pi * t) - par.S2 * x^2)
solar(x::Vec, t::Float64, ::Val{:water}, par::Par)::Vec =
    @. (par.a0 - par.a2 * x^2) * (par.S0 - par.S1 * x * cos(2pi * t) - par.S2 * x^2)
solar(x::Vec, t::Float64, surface::Symbol, par::Par)::Vec = solar(x, t, Val(surface), par)

# phi-weighted average
weighted_avg(vi::Vec, vw::Vec, phi::Vec)::Vec = @. vi*phi + (1-phi)vw

# temperatures
water_temp(Ew::Vec, par::Par)::Vec = @. par.Tm + Ew / par.cw

ice_temp(T0::Vec, par::Par)::Vec = min.(T0, par.Tm)

solveT0(x::Vec, t::Float64, h::Vec, Tg::Vec, Tw::Vec, phi::Vec, f::Float64, par::Par)::Vec =
    @. (
        $(solar(x, t, :ice, par)) - par.A + f - (1-phi)Tw * (par.B + par.cg/par.tau)
        + par.Tm * (par.B + par.k/h) + par.cg/par.tau * Tg
    ) / (phi * (par.B + par.cg/par.tau) + par.k/h)

function stepTg!(
    t::Float64, Tg::Vec, h::Vec, T0::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime, par::Par
)::Vec
    frez = @. (T0<par.Tm) & (h>0)
    watr = .~frez
    diagphi = SA.spdiagm(phi)
    invM = SA.spdiagm(inv.(phi * (par.B + par.cg/par.tau) .+ par.k./h))
    Tg .= (
            (1+st.dt/par.tau)LA.I
            - st.dt*par.D/par.cg * get_diffop(st)
            - (st.dt*par.cg/par.tau^2 * diagphi * invM)SA.spdiagm(frez)
        ) \ (
            Tg
            + st.dt/par.tau * (
                (LA.I-diagphi)Tw
                + (
                    diagphi * (
                        par.Tm * (par.B*LA.I + par.k*LA.I * SA.spdiagm(inv.(h)))
                        - (par.B + par.cg/par.tau) * (LA.I-diagphi)SA.spdiagm(Tw)
                        + SA.spdiagm(solar(st.x, t, :ice, par)) - par.A*LA.I + f*LA.I
                    ) * invM
                )frez
                + par.Tm * diagphi*watr
            )
        )
    return Tg
end # function stepTg!

# lateral melt rate
wlat(Tw::Vec, par::Par)::Vec = @. par.m1 * (Tw - par.Tm)^par.m2

# concentration
function concentration(Ei::Vec, h::Vec, par::Par)::Vec
    phi = @. -Ei / (par.Lf * h)
    zeroref!(phi, h)
    condset!(phi, 1.0, >(1)) # correct concentration
    # phi = @. Float64(Ei<0) # reproducing WE15
end # function concentration

# floe number
function num(D::Vec, phi::Vec, par::Par)::Vec
    n = @. phi / (par.alpha * D^2)
    zeroref!(n, D)
    return n
end # function num

# lead region area
function area_lead(D::Vec, phi::Vec, n::Vec, par::Par)::Vec
    ring = @. par.alpha * n * ((D + 2par.rl)^2 - D^2)
    return min.(ring, 1 .- phi)
end # function area_lead

# fluxes
function vert_flux(
    t::Float64, surface::Symbol, Tg::Vec, Tbar::Vec, f::Float64, st::SpaceTime, par::Par
)::Vec
    L = @. par.A + par.B * (Tbar - par.Tm) # OLR
    return solar(st.x, t, surface, par) .- L .+ par.cg/par.tau * (Tg-Tbar) .+ par.Fb .+ f
end # function vert_flux

function lat_flux(h::Vec, D::Vec, Tw::Vec, phi::Vec, par::Par)::Vec
    Flat = @. phi * h * par.Lf * $(wlat(Tw, par)) * pi / (par.alpha*D)
    zeroref!(Flat, D)
    return Flat
end # function lat_flux

bot_flux(Tw::Vector, par::Collection) = par.rhow * par.cp * par.ch * par.u0 * (Tw .- par.Tm)

function redistributeE(rEi::Vec, rEw::Vec)::NTuple{4,Vec}
    cEi = clamp.(rEi, -Inf, 0)
    cEw = clamp.(rEw, 0, Inf)
    psiEidt = rEi .- cEi # +
    psiEwdt = rEw .- cEw # -
    Ei = cEi .+ psiEwdt # -
    Ew = cEw .+ psiEidt # +
    return (Ei, Ew, psiEidt, psiEwdt)
end # function redistributeE

# redistribution functions
function split_psiEw(psiEw::Vec, phi::Vec, Al::Vec)::Tuple{Vec,Vec}
    Ql = @. Al / (1-phi) * psiEw
    condset!(Ql, 0.0, isone, phi) # fix rounding errors
    Qp = psiEw - Ql
    return (Ql, Qp)
end # function split_psiEw

dphip(Qp::Vec, par::Par)::Vec = @. -Qp / (par.Lf * par.hmin) # change rate of φ due to pancakes

function average(f::Vec, fn::Float64, n::Vec, dn::Vec)::Vec
    total = n .+ dn
    avgd = @. (n*f + dn*fn) / total
    zeroref!(avgd, total)
    # return f # reproducing WE15
    return avgd
end # function average

# differential equations
Ei_t(phi::Vec, Fvi::Vec, Flat::Vec, Fbot::Vec)::Vec = @. phi * Fvi + Flat + phi * Fbot
Ew_t(phi::Vec, Fvw::Vec, Flat::Vec, Fbot::Vec)::Vec = @. (1-phi)Fvw - Flat - phi * Fbot
h_t(Fvi::Vec, Fbot::Vec, par::Par)::Vec = -1/par.Lf * (Fvi + Fbot)
function D_t(h::Vec, D::Vec, Ti::Vec, Tw::Vec, phi::Vec, Ql::Vec, par::Par; breakup::BitArray)::Vec
    lat_melt = -pi / 2par.alpha * wlat(Tw, par)
    lat_grow = @. -D / (2 * par.Lf * h * phi) * Ql
    weld = @. par.kappa * par.alpha / 4 * phi * D^3
    zeroref!(lat_grow, h)
    condset!(weld, 0.0, >=(par.Tm), Ti)
    weld[breakup] .= 0.0 # no welding if breaking
    return @. lat_melt + lat_grow + weld
end # function D_t

forward_euler(var::Vec, grad::Vec, dt::Float64)::Vec = @. var + grad*dt

# common template of initialise function for MIZModel and WIModel
function _initialise(
    model::Union{MIZModel,WIModel}, st::SpaceTime, forcing::Forcing, par::Par, init::Collection{Vec};
    lastonly::Bool
) # -> Tuple{Collection{Vec}, Solutions{M,F,V}, Solutions{M,F,V}}
    # create storages
    solvars = Set{Symbol}((:Ei, :Ew, :D, :h, :E, :Ti, :Tw, :T, :phi, :n))
    vars, sols, annusol = create_storages(model, solvars, st, forcing, par, init; lastonly)
    # diagnostic variables read by the first step!, on the same timestep as Ei, Ew, h and D
    vars.phi = concentration(vars.Ei, vars.h, par)
    vars.Tw = water_temp(vars.Ew, par)
    vars.T0 = solveT0(st.x, st.T[1], vars.h, vars.Tg, vars.Tw, vars.phi, forcing(st.T[1]), par)
    vars.Ti = condset(ice_temp(vars.T0, par), 0.0, isnan) # eliminate NaNs for calculations
    vars.T = weighted_avg(vars.Ti, vars.Tw, vars.phi)
    return (vars, sols, annusol)
end # function _initialise

Infrastructure.initialise(
    model::MIZModel, st::SpaceTime, forcing::Forcing, par::Par, init::Collection{Vec};
    lastonly::Bool=true
) = _initialise(model, st, forcing, par, init; lastonly)
    # -> Tuple{Collection{Vec}, Solutions{MIZModel,F,V}, Solutions{MIZModel,F,V}}

function Infrastructure.step!(
    ::MIZModel, t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime, par::Par;
    breakup::BitArray=falses(st.nx)
)::Collection{Vec}
    # prepare for the step
    vars.n = num(vars.D, vars.phi, par) # refresh in case D was changed by wave breakup
    Ti = condset(vars.Ti, 0.0, isnan) # eliminate NaNs for calculations
    # calculate fluxes
    Fvi = vert_flux(t, :ice, vars.Tg, vars.T, f, st, par)
    Fvw = vert_flux(t, :water, vars.Tg, vars.T, f, st, par)
    Flat = lat_flux(vars.h, vars.D, vars.Tw, vars.phi, par)
    Fbot = bot_flux(vars.Tw, par)
    # update enthalpy
    rEi = forward_euler(vars.Ei, Ei_t(vars.phi, Fvi, Flat, Fbot), st.dt)
    rEw = forward_euler(vars.Ew, Ew_t(vars.phi, Fvw, Flat, Fbot), st.dt)
    Ei, Ew, _, psiEwdt = redistributeE(rEi, rEw)
    vars.Ei = Ei # !
    vars.Ew = Ew # !
    vars.E = vars.Ei + vars.Ew # !
    # update floe size and thickness
    Al = area_lead(vars.D, vars.phi, vars.n, par)
    Ql, Qp = split_psiEw(psiEwdt/st.dt, vars.phi, Al)
    phip = st.dt * dphip(Qp, par)
    lasth = vars.h # save for D
    vars.h = forward_euler(
        average(vars.h, par.hmin, vars.phi, phip), # new pancakes
        h_t(Fvi, Fbot, par),
        st.dt
    ) # !
    vars.D = forward_euler(
        average(vars.D, par.Dmin, vars.phi, phip), # new pancakes
        D_t(lasth, vars.D, Ti, vars.Tw, vars.phi, Ql, par; breakup),
        st.dt
    ) # !
    clamp!(vars.h, 0, Inf) # avoid overshooting to negative thickness
    zeroref!(vars.h, vars.Ei) # restrict non-existence
    clamp!(vars.D, 0, par.Dmax)
    zeroref!(vars.D, vars.Ei) # restrict non-existence
    # advance the diagnostic variables
    vars.phi = concentration(vars.Ei, vars.h, par) # !
    vars.n = num(vars.D, vars.phi, par) # !
    vars.Tw = water_temp(vars.Ew, par) # !
    vars.Tg = stepTg!(t+st.dt, vars.Tg, vars.h, vars.T0, vars.Tw, vars.phi, f, st, par) # !
    vars.T0 = solveT0(st.x, t+st.dt, vars.h, vars.Tg, vars.Tw, vars.phi, f, par)
    vars.Ti = condset(ice_temp(vars.T0, par), 0.0, isnan) # !
    vars.T = weighted_avg(vars.Ti, vars.Tw, vars.phi) # !
    # set NaNs to no existence
    condset!(vars.Ti, NaN, iszero, vars.Ei)
    return vars
end # function Infrastructure.step!

end # module MIZEBM
