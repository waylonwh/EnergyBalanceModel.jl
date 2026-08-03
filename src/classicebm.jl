module ClassicEBM # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import LinearAlgebra as LA, SparseArrays as SA
import NonlinearSolve as NlinSol

@persistent(
    cg_tau::Float64, dt_tau::Float64, dc::Float64, kappa::Matrix{Float64},
    S::Matrix{Float64}, M::Float64, aw::Vec, kLf::Float64,
    id::UInt = UInt(0),
    @inline function get_statics(st::SpaceTime, par::Collection)::@NamedTuple{
        cg_tau::Float64, dt_tau::Float64, dc::Float64, kappa::Matrix{Float64},
        S::Matrix{Float64}, M::Float64, aw::Vec, kLf::Float64
    }
        if id != hash((st, par)) # recompute only if st or par changed
            # Difinitions for implicit scheme for Tg
            cg_tau = par.cg / par.tau
            dt_tau = st.dt / par.tau
            dc = dt_tau * cg_tau
            kappa = (1+dt_tau) * LA.I(st.nx) - st.dt * par.D * get_diffop(st) / par.cg
            # Seasonal forcing [WE15 Eq. (3)]
            S = repeat(par.S0 .- par.S2 * st.x.^2, 1, st.nt) -
                repeat(par.S1 * cos.(2pi*st.t'), st.nx, 1) .* repeat(st.x, 1, st.nt)
            S = hcat(S, S[:,1])
            # Further definitions
            M = par.B + cg_tau
            aw = @. par.a0 - par.a2 * st.x^2
            kLf = par.k * par.Lf
            # update id
            id = hash((st, par))
        end # if !=
        return (; cg_tau, dt_tau, dc, kappa, S, M, aw, kLf)
    end # function get_statics
) # @persistent

surface_temperature(E::Vector, T0::Vector, par::Collection) = # -> Vector
    @. E/par.cw * (E>=0) + T0 * (E<0)*(T0<0) # WE15 Eq. (9)

T0eq(T0::Vector, h::Vector, aS::Vector, diffop::AbstractMatrix, f::Real, par::Collection) = # -> Vector
    @. par.k * (par.Tm-T0) / h + aS - par.A - par.B * (T0-par.Tm) + par.D * diffop*T0 + f

function solveT0(
    guess::Vector, h::Vector, aS::Vector, diffop::AbstractMatrix, f::Real, par::Collection;
    abstol::AbstractFloat=1e-10
) # -> Vector
    prob = NlinSol.NonlinearProblem{false}(
        (T0, p) -> T0eq(T0, p.h, p.aS, p.diffop, p.f, p.par),
        guess, (; h, aS, diffop, f, par)
    )
    sol = NlinSol.solve(prob, NlinSol.TrustRegion(); abstol)
    NlinSol.SciMLBase.successful_retcode(sol) ||
        @warn(
            "Nonlinear solver did not converge when solving for wavenumber. Result may be inaccurate.",
            sol.retcode, sol.resid, sol.u
        )
    return sol.u
end # function solveT0

Infrastructure.initialise(
    model::ClassicModel, st::SpaceTime, forcing::Forcing, par::Collection, init::Collection{Vec};
    lastonly::Bool=true
) = create_storages(model, Set{Symbol}((:E, :T, :h)), st, forcing, par, init; lastonly)
    # -> Tuple{Collection{Vec},Solutions{ClassicModel,F,V},Solutions{ClassicModel,F,V}}

function specialised_step!(
    ::GhostLayerSolver, vars::Collection{<:Vector}, C::Vector, stat::NamedTuple,
    st::SpaceTime, par::Collection, _...
)::Nothing
    # surface temperature
    T0 = @. C / (stat.M - stat.kLf/vars.E) # WE15 Eq. (A3)
    vars.T = surface_temperature(vars.E, T0, par)
    # Forward Euler for E
    @. vars.E += st.dt * (C - stat.M*vars.T + par.Fb) # WE15 Eq. (A2)
    # Implicit Euler for Tg
    vars.Tg =
        (stat.kappa - SA.spdiagm(stat.dc ./ (stat.M .- stat.kLf./vars.E) .* (T0.<0).*(vars.E.<0))) \
        (
            vars.Tg +
            (
                stat.dt_tau * (vars.E/par.cw.*(vars.E.>=0) +
                (par.ai*view(stat.S, :, i+1) .- par.A .+ f) ./ (stat.M .- stat.kLf./vars.E) .* (T0.<0).*(vars.E.<0))
            )
        ) # \ # WE15 Eq. (A1)
    return nothing
end # function specialised_step!

function specialised_step!(
    ::NonlinearSolver, vars::Collection{Vector{FT}}, C::Vector{FT}, stat::NamedTuple,
    st::SpaceTime, par::Collection{FT}, aS::Vector{FT}, f::FT
)::Nothing where FT<:AbstractFloat
    vars.T0 = solveT0(
        get(vars, :T0, zeros(FT, length(C))), vars.h, aS, get_diffop(st), f, par
    )
    vars.T = surface_temperature(vars.E, vars.T0, par)
    @. vars.E += st.dt * (C - stat.M*vars.T + par.Fb) # WE15 Eq. (A2)
    return nothing
end # function specialised_step!

function specialised_step!(::ActiveSetSolver, vars::Collection)::Nothing
    return nothing
end # function specialised_step!

function Infrastructure.step!(
    ::ClassicModel, t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime, par::Collection;
    solver::AbstractSolver=GhostLayerSolver()
)::Collection{Vec}
    # get static variables
    stat = get_statics(st, par)
    # get time index
    i = round(Int, mod1((t + st.dt/2) * st.nt, st.nt))
    # forcing
    alpha = @. stat.aw * (vars.E>=0) + par.ai * (vars.E<0) # WE15 Eq. (4)
    C = @. alpha*stat.S[:,i] + stat.cg_tau*vars.Tg - par.A + f
    # step T and E
    specialised_step!(solver, vars, C, stat, st, par, alpha*stat.S[:,i], f)
    # Infer ice thickness
    vars.h = @. -vars.E / par.Lf * (vars.E<0)
    return vars
end # function Infrastructure.step!

end # module ClassicEBM
