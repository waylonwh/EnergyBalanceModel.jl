module MIZEBM # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import LinearAlgebra as LA
import SparseArrays as SA

# solar radiation absorbed on ice and water
function solar(x::Vector{T}, t::T, ::Val{:ice}, par::Par{T})::Vector{T} where {T<:AbstractFloat}
    return @. par.ai * (par.S0 - par.S1 * x * cos(T(2pi) * t) - par.S2 * x^2)
end

function solar(x::Vector{T}, t::T, ::Val{:water}, par::Par{T})::Vector{T} where {T<:AbstractFloat}
    return @. (par.a0 - par.a2 * x^2) * (par.S0 - par.S1 * x * cos(T(2pi) * t) - par.S2 * x^2)
end

function solar(x::Vector{T}, t::T, surface::Symbol, par::Par{T})::Vector{T} where {T<:AbstractFloat}
    return solar(x, t, Val(surface), par)
end

# phi-weighted average
function weighted_avg(vi::Vector{T}, vw::Vector{T}, phi::Vector{T}) where {T<:AbstractFloat}
    return @. vi*phi + (1-phi)vw
end

# temperatures
function water_temp(Ew::Vector{T}, phi::Vector{T}, par::Par{T}) where {T<:AbstractFloat}
    return @. par.Tm + Ew / ((1-phi)par.cw)
end

function water_temp_nonan(Ew::Vector{T}, phi::Vector{T}, par::Par{T}) where {T<:AbstractFloat}
    return condset!(water_temp(Ew, phi, par), zero(T), isone, phi)
end

function ice_temp(T0::Vector{T}, par::Par{T}) where {T<:AbstractFloat}
    return min.(T0, par.Tm)
end

function solveT0(
    x::Vector{T}, t::T, h::Vector{T}, Tg::Vector{T}, Tw::Vector{T}, phi::Vector{T}, f::T, par::Par{T}
) where {T<:AbstractFloat}
    return @. (
        $(solar(x, t, :ice, par)) - par.A + f - (1-phi)Tw * (par.B + par.cg/par.tau)
        + par.Tm * (par.B + par.k/h) + par.cg/par.tau * Tg
    ) / (phi * (par.B + par.cg/par.tau) + par.k/h)
end

function stepTg!(
    t::T, Tg::Vector{T}, h::Vector{T}, T0::Vector{T}, Tw::Vector{T}, phi::Vector{T}, f::T, st::SpaceTime{F,T}, par::Par{T}
)::Vector{T} where {F,T<:AbstractFloat}
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
function wlat(Tw::Vector{T}, par::Par{T}) where {T<:AbstractFloat}
    return @. par.m1 * (Tw - par.Tm^par.m2)
end

# concentration
function concentration(Ei::Vector{T}, h::Vector{T}, par::Par{T})::Vector{T} where {T<:AbstractFloat}
    phi = @. -Ei / (par.Lf * h)
    zeroref!(phi, h)
    condset!(phi, one(T), >(one(T))) # correct concentration
    # phi = @. Float64(Ei<0) # reproducing WE15
end # function concentration

# floe number
function num(D::Vector{T}, phi::Vector{T}, par::Par{T})::Vector{T} where {T<:AbstractFloat}
    n = @. phi / (par.alpha * D^2)
    zeroref!(n, D)
    return n
end # function num

# lead region area
function area_lead(D::Vector{T}, phi::Vector{T}, n::Vector{T}, par::Par{T})::Vector{T} where {T<:AbstractFloat}
    ring = @. par.alpha * n * ((D + 2par.rl)^2 - D^2)
    return min.(ring, one(T) .- phi)
end # function area_lead

# fluxes
function vert_flux(
    t::T, surface::Symbol, Tg::Vector{T}, Tbar::Vector{T}, f::T, st::SpaceTime{F,T}, par::Par{T}
)::Vector{T} where {F,T<:AbstractFloat}
    L = @. par.A + par.B * (Tbar - par.Tm) # OLR
    return solar(st.x, t, surface, par) .- L .+ par.cg/par.tau * (Tg-Tbar) .+ par.Fb .+ f
end # function vert_flux

function lat_flux(h::Vector{T}, D::Vector{T}, Tw::Vector{T}, phi::Vector{T}, par::Par{T})::Vector{T} where {T<:AbstractFloat}
    Flat = @. phi * h * par.Lf * $(wlat(Tw, par)) * T(pi) / (par.alpha*D)
    zeroref!(Flat, D)
    return Flat
end # function lat_flux

function redistributeE(rEi::Vector{T}, rEw::Vector{T})::NTuple{4,Vector{T}} where {T<:AbstractFloat}
    cEi = clamp.(rEi, -T(Inf), zero(T))
    cEw = clamp.(rEw, zero(T), T(Inf))
    psiEidt = rEi .- cEi # +
    psiEwdt = rEw .- cEw # -
    Ei = cEi .+ psiEwdt # -
    Ew = cEw .+ psiEidt # +
    return (Ei, Ew, psiEidt, psiEwdt)
end # function redistributeE

# redistribution functions
function split_psiEw(psiEw::Vector{T}, phi::Vector{T}, Al::Vector{T})::Tuple{Vector{T},Vector{T}} where {T<:AbstractFloat}
    Ql = @. Al / (1-phi) * psiEw
    condset!(Ql, zero(T), isone, phi) # fix rounding errors
    Qp = psiEw - Ql
    return (Ql, Qp)
end # function split_psiEw

function dphip(Qp::Vector{T}, par::Par{T}) where {T<:AbstractFloat}
    return @. -Qp / (par.Lf * par.hmin) # change rate of φ due to pancakes
end

function average(f::Vector{T}, fn::T, n::Vector{T}, dn::Vector{T})::Vector{T} where {T<:AbstractFloat}
    total = n .+ dn
    avgd = @. (n*f + dn*fn) / total
    zeroref!(avgd, total)
    # return f # reproducing WE15
    return avgd
end # function average

# differential equations
function Ei_t(phi::Vector{T}, Fvi::Vector{T}, Flat::Vector{T}) where {T<:AbstractFloat}
    return @. phi * Fvi + Flat
end

function Ew_t(phi::Vector{T}, Fvw::Vector{T}, Flat::Vector{T}) where {T<:AbstractFloat}
    return @. (1-phi)Fvw - Flat
end

function h_t(Fvi::Vector{T}, par::Par{T}) where {T<:AbstractFloat}
    return -one(T)/par.Lf * Fvi
end
function D_t(h::Vector{T}, D::Vector{T}, Ti::Vector{T}, Tw::Vector{T}, phi::Vector{T}, Ql::Vector{T}, par::Par{T})::Vector{T} where {T<:AbstractFloat}
    lat_melt = -T(pi) / T(2) * par.alpha * wlat(Tw, par)
    lat_grow = @. -D / (2 * par.Lf * h * phi) * Ql
    weld = @. par.kappa * par.alpha / 4 * phi * D^3
    zeroref!(lat_grow, h)
    condset!(weld, zero(T), >=(par.Tm), Ti)
    return @. lat_melt + lat_grow + weld
end # function D_t

function forward_euler(var::Vector{T}, grad::Vector{T}, dt::T) where {T<:AbstractFloat}
    return @. var + grad*dt
end

# common template of initialise function for MIZModel and WIModel
function _initialise(
    model::Union{MIZModel,WIModel},
    st::SpaceTime{F,T},
    forcing::Forcing{T,C},
    par::Par{T},
    init::Collection{Vector{T}};
    lastonly::Bool
) where {F,C,T<:AbstractFloat} # -> Tuple{Collection{Vector}, Solutions{M,F,V}, Solutions{M,F,V}}
    # create storages
    solvars = Set{Symbol}((:Ei, :Ew, :D, :h, :E, :Ti, :Tw, :T, :phi, :n))
    model isa WIModel && push!(solvars, :Ewave, :lambda) # add wave variables for WIModel
    vars, sols, annusol = create_storages(model, solvars, st, forcing, par, init; lastonly)
    # compute phi and Tw
    vars.nextphi = concentration(vars.Ei, vars.h, par)
    vars.nextTw = water_temp(vars.Ew, vars.nextphi, par)
    vars.nextT0 = solveT0(st.x, st.T[1], vars.h, vars.Tg, vars.nextTw, vars.nextphi, forcing(st.T[1]), par)
    condset!(vars.nextTw, zero(T), isnan) # eliminate NaNs for calculations
    return (vars, sols, annusol)
end # function _initialise

Infrastructure.initialise(
    model::MIZModel, st::SpaceTime{F,T}, forcing::Forcing{T,C}, par::Par{T}, init::Collection{Vector{T}};
    lastonly::Bool=true
) where {F,C,T<:AbstractFloat} = _initialise(model, st, forcing, par, init; lastonly)
    # -> Tuple{Collection{Vector}, Solutions{MIZModel,F,V}, Solutions{MIZModel,F,V}}

function Infrastructure.step!(
    ::MIZModel, t::T, f::T, vars::Collection{Vector{T}}, st::SpaceTime{F,T}, par::Par{T}
)::Collection{Vector{T}} where {F,T<:AbstractFloat}
    # assign next variables to current
    vars.phi = vars.nextphi
    vars.Tw = vars.nextTw
    T0 = vars.nextT0
    # compute diagnostic variables
    vars.Ti = ice_temp(T0, par)
    condset!(vars.Ti, zero(T), isnan) # eliminate NaNs for calculations
    vars.T = weighted_avg(vars.Ti, vars.Tw, vars.phi)
    vars.n = num(vars.D, vars.phi, par)
    # calculate fluxes
    Fvi = vert_flux(t, :ice, vars.Tg, vars.T, f, st, par)
    Fvw = vert_flux(t, :water, vars.Tg, vars.T, f, st, par)
    Flat = lat_flux(vars.h, vars.D, vars.Tw, vars.phi, par)
    # update enthalpy
    rEi = forward_euler(vars.Ei, Ei_t(vars.phi, Fvi, Flat), st.dt)
    rEw = forward_euler(vars.Ew, Ew_t(vars.phi, Fvw, Flat), st.dt)
    Ei, Ew, _, psiEwdt = redistributeE(rEi, rEw)
    vars.Ei = Ei # !
    vars.Ew = Ew # !
    vars.E = weighted_avg(vars.Ei, vars.Ew, vars.phi) # !
    # update floe size and thickness
    Al = area_lead(vars.D, vars.phi, vars.n, par)
    Ql, Qp = split_psiEw(psiEwdt/st.dt, vars.phi, Al)
    phip = st.dt * dphip(Qp, par)
    lasth = vars.h # save for D
    vars.h = forward_euler(
        average(vars.h, par.hmin, vars.phi, phip), # new pancakes
        h_t(Fvi, par),
        st.dt
    ) # !
    vars.D = forward_euler(
        average(vars.D, par.Dmin, vars.phi, phip), # new pancakes
        D_t(lasth, vars.D, vars.Ti, vars.Tw, vars.phi, Ql, par),
        st.dt
    ) # !
    clamp!(vars.h, zero(T), T(Inf)) # avoid overshooting to negative thickness
    zeroref!(vars.h, vars.Ei) # restrict non-existence
    clamp!(vars.D, zero(T), par.Dmax)
    zeroref!(vars.D, vars.Ei) # restrict non-existence
    # update variables for Tg
    vars.nextphi = concentration(vars.Ei, vars.h, par) # !
    vars.nextTw = water_temp_nonan(vars.Ew, vars.nextphi, par) # !
    vars.nextT0 = solveT0(st.x, t, vars.h, vars.Tg, vars.nextTw, vars.nextphi, f, par)
    vars.Tg = stepTg!(t, vars.Tg, vars.h, vars.nextT0, vars.nextTw, vars.nextphi, f, st, par) # !
    # set NaNs to no existence
    condset!(vars.Ti, T(NaN), iszero, vars.Ei)
    condset!(vars.Tw, T(NaN), >(T(0.95)), vars.phi)
    return vars
end # function Infrastructure.step!

end # module MIZEBM
