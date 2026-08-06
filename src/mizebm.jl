module MIZEBM # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import LinearAlgebra as LA
import SparseArrays as SA
import NonlinearSolve as NlinSol

# solar radiation absorbed on ice and water
solar(x::Vec, t::Float64, ::Val{:ice}, par::Collection)::Vec =
    @. par.ai * (par.S0 - par.S1 * x * cos(2pi * t) - par.S2 * x^2)
solar(x::Vec, t::Float64, ::Val{:water}, par::Collection)::Vec =
    @. (par.a0 - par.a2 * x^2) * (par.S0 - par.S1 * x * cos(2pi * t) - par.S2 * x^2)
solar(x::Vec, t::Float64, surface::Symbol, par::Collection)::Vec = solar(x, t, Val(surface), par)

# phi-weighted average
weighted_avg(vi::Vector, vw::Vector, phi::Vector) = @. vi*phi + (1-phi)vw # -> Vector

# temperatures
water_temp(Ew::Vec, par::Collection)::Vec = @. par.Tm + Ew / par.cw
ice_temp(T0::Vector, par::Collection) = min.(T0, par.Tm) # -> Vector

# ghost layer scheme
solveT0(
    ::GhostLayerSolver, x::Vector, t::Float64, h::Vector, Tg::Vector, Tw::Vector,
    phi::Vector, f::Real, par::Collection
) = @. ( # -> Vector
    $(solar(x, t, :ice, par)) - par.A + f - (1-phi)Tw * (par.B + par.cg/par.tau)
    + par.Tm * (par.B + par.k/h) + par.cg/par.tau * Tg
) / (phi * (par.B + par.cg/par.tau) + par.k/h)

function stepTg!(
    t::Float64, Tg::Vec, h::Vec, T0::Vec, Tw::Vec, phi::Vec, f::Float64, st::SpaceTime, par::Collection
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

# nonlinear scheme
function T0eq(
    T0::Vector, x::Vector, t::Real, h::Vector, Tw::Vector, phi::Vector, f::Real,
    diffop::AbstractMatrix, par::Collection
) # -> Vector
    T = weighted_avg(ice_temp(T0, par), Tw, phi)
    vec = @. par.k * (par.Tm - T0) / h # SCM
    vec .+= solar(x, t, Val(:ice), par) # solar on ice
    @. vec -= par.A + par.B * (T - par.Tm) # OLR
    vec .+= par.D * diffop * T # diffusion
    vec .+= f # forcing
    return vec
end # function T0eq

function solveT0(
    ::NonlinearSolver, guess::Vec, x::Vector, t::Real, h::Vector, Tw::Vector, phi::Vector,
    diffop::AbstractMatrix, f::Real, par::Collection; abstol::AbstractFloat=1e-10
)::Vec
    hp = replace(h, 0.0=>par.hmin) # avoid division by zero where there is no ice
    init = replace(guess, NaN=>par.Tm) # NaNs cannot be used as an initial guess
    prob = NlinSol.NonlinearProblem{false}(
        (T0, p) -> T0eq(T0, p.x, p.t, p.hp, p.Tw, p.phi, p.f, p.diffop, p.par),
        init, (; x, t, hp, Tw, phi, f, diffop, par)
    )
    sol = NlinSol.solve(prob, NlinSol.NewtonRaphson(); abstol)
    NlinSol.SciMLBase.successful_retcode(sol) ||
        @warn(
            "Nonlinear solver did not converge when solving for surface temperature. Result may be inaccurate.",
            sol.retcode, sol.resid, sol.u
        )
    return sol.u
end # function solveT0

# lateral melt rate
wlat(Tw::Vec, par::Collection)::Vec = @. par.m1 * (Tw - par.Tm)^par.m2

# concentration
function concentration(Ei::Vec, h::Vec, par::Collection)::Vec
    phi = @. -Ei / (par.Lf * h)
    zeroref!(phi, h)
    replace!(p -> p>1.0 ? 1.0 : p, phi) # correct concentration
    # phi = @. Float64(Ei<0) # reproducing WE15
end # function concentration

# floe number
function num(D::Vec, phi::Vec, par::Collection)::Vec
    n = @. phi / (par.alpha * D^2)
    zeroref!(n, D)
    return n
end # function num

# lead region area
function area_lead(D::Vec, phi::Vec, n::Vec, par::Collection)::Vec
    ring = @. par.alpha * n * ((D + 2par.rl)^2 - D^2)
    return min.(ring, 1 .- phi)
end # function area_lead

# fluxes
function vert_flux(
    ::GhostLayerSolver, t::Float64, surface::Symbol, Tbar::Vector, f::Real, st::SpaceTime,
    par::Collection, Tg::Vector
) # -> Vector
    L = @. par.A + par.B * (Tbar - par.Tm) # OLR
    return solar(st.x, t, surface, par) .- L .+ par.cg/par.tau * (Tg-Tbar) .+ par.Fb .+ f
end # function vert_flux

function vert_flux(
    ::NonlinearSolver, t::Float64, surface::Symbol, Tbar::Vector, f::Real, st::SpaceTime,
    par::Collection, _
) # -> Vector
    L = @. par.A + par.B * (Tbar - par.Tm) # OLR
    return solar(st.x, t, surface, par) .- L .+ par.D * get_diffop(st) * Tbar .+ par.Fb .+ f
end # function vert_flux

function lat_flux(h::Vec, D::Vec, Tw::Vec, phi::Vec, par::Collection)::Vec
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

dphip(Qp::Vec, par::Collection)::Vec = @. -Qp / (par.Lf * par.hmin) # change rate of φ due to pancakes

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
h_t(Fvi::Vec, Fbot::Vec, par::Collection)::Vec = -1/par.Lf * (Fvi + Fbot)
function D_t(h::Vec, D::Vec, Ti::Vec, Tw::Vec, phi::Vec, Ql::Vec, par::Collection; breakup::BitArray)::Vec
    lat_melt = -pi / 2par.alpha * wlat(Tw, par)
    lat_grow = @. -D / (2 * par.Lf * h * phi) * Ql
    weld = @. par.kappa * par.alpha / 4 * phi * D^3
    zeroref!(lat_grow, h)
    condset!(weld, 0.0, >=(par.Tm), Ti)
    weld[breakup] .= 0.0 # no welding if breaking
    return @. lat_melt + lat_grow + weld
end # function D_t

forward_euler(var::Vec, grad::Vec, dt::Float64)::Vec = @. var + grad*dt

function step_temperature!(
    solver::GhostLayerSolver, vars::Collection{<:Vector}, t::Real, st::SpaceTime,
    par::Collection, f::Real; initstep::Bool=false
) # -> Vector
    initstep || (vars.Tg = stepTg!(t, vars.Tg, vars.h, vars.T0, vars.Tw, vars.phi, f, st, par))
    vars.T0 = solveT0(solver, st.x, t, vars.h, vars.Tg, vars.Tw, vars.phi, f, par)
    return vars.T0
end # function specialised_step!

step_temperature!(
    solver::NonlinearSolver, vars::Collection{<:Vector}, t::Real, st::SpaceTime,
    par::Collection, f::Real; _...
) = vars.T0=solveT0(
    solver, get(vars, :T0, fill(par.Tm, st.nx)),
    st.x, t, vars.h, vars.Tw, vars.phi, get_diffop(st), f, par
) # -> Vector

function step_temperature!(::ActiveSetSolver, vars::Collection)
    throw(ErrorException("ActiveSetSolver is not implemented for MIZModel or WIModel."))
end # function specialised_step!

# common template of initialise function for MIZModel and WIModel
function _initialise(
    model::Union{MIZModel,WIModel}, st::SpaceTime, forcing::Forcing, par::Collection, init::Collection{Vec};
    solver::AbstractSolver, lastonly::Bool
) # -> Tuple{Collection{Vec}, Solutions{M,F,V}, Solutions{M,F,V}}
    # create storages
    solvars = Set{Symbol}((:Ei, :Ew, :D, :h, :E, :Ti, :Tw, :T, :phi, :n))
    vars, sols, annusol = create_storages(model, solvars, st, forcing, par, init; solver, lastonly)
    # diagnostic variables read by the first step!, on the same timestep as Ei, Ew, h and D
    vars.phi = concentration(vars.Ei, vars.h, par)
    vars.Tw = water_temp(vars.Ew, par)
    step_temperature!(solver, vars, st.T[1], st, par, forcing(st.T[1]); initstep=true) # step T0
    vars.Ti = replace!(ice_temp(vars.T0, par), NaN=>0.0) # eliminate NaNs for calculations
    vars.T = weighted_avg(vars.Ti, vars.Tw, vars.phi)
    return (vars, sols, annusol)
end # function _initialise

Infrastructure.initialise(
    model::MIZModel, st::SpaceTime, forcing::Forcing, par::Collection, init::Collection{Vec};
    solver::AbstractSolver, lastonly::Bool, _...
) = _initialise(model, st, forcing, par, init; solver, lastonly)
    # -> Tuple{Collection{Vec}, Solutions{MIZModel,F,V}, Solutions{MIZModel,F,V}}

function Infrastructure.step!(
    ::MIZModel, t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime, par::Collection;
    solver::AbstractSolver, breakup::BitArray=falses(st.nx), _...
)::Collection{Vec}
    # prepare for the step
    vars.n = num(vars.D, vars.phi, par) # refresh in case D was changed by wave breakup
    Ti = replace(vars.Ti, NaN=>0.0) # eliminate NaNs for calculations
    # calculate fluxes
    Fvi = vert_flux(solver, t, :ice, vars.T, f, st, par, vars.Tg)
    Fvw = vert_flux(solver, t, :water, vars.T, f, st, par, vars.Tg)
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
    # step the surface temperature; time step n+1 start point
    step_temperature!(solver, vars, t+st.dt, st, par, f) # ! step T0
    vars.Ti = ice_temp(vars.T0, par) # !
    vars.T = weighted_avg(replace(vars.Ti, NaN=>0.0), vars.Tw, vars.phi)
    # set NaNs to no existence
    condset!(vars.Ti, NaN, iszero, vars.Ei)
    return vars
end # function Infrastructure.step!

end # module MIZEBM
