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

T0eq(T0::Vector, E::Vector, h::Vector, aS::Vector, diffop::AbstractMatrix, f::Real, par::Collection) = # -> Vector
    @. ifelse(
        E < 0,
        par.k * (0-T0) / h + aS - par.A - par.B * (T0-0) + par.D * $(diffop*surface_temperature(E, T0, par)) + f,
        T0 - 0
    )

function solveT0(
    guess::Vector, E::Vector, h::Vector, aS::Vector, diffop::AbstractMatrix, f::Real, par::Collection;
    abstol::AbstractFloat=1e-10
) # -> Vector
    prob = NlinSol.NonlinearProblem{false}(
        (T0, p) -> T0eq(T0, p.E, p.h, p.aS, p.diffop, p.f, p.par),
        guess, (; E, h, aS, diffop, f, par)
    )
    sol = NlinSol.solve(prob, NlinSol.NewtonRaphson(); abstol)
    NlinSol.SciMLBase.successful_retcode(sol) ||
        @warn(
            "Nonlinear solver did not converge when solving for wavenumber. Result may be inaccurate.",
            sol.retcode, sol.resid, sol.u
        )
    return sol.u
end # function solveT0

function solveT(
    vars::Collection{<:Vector}, aS::Vector, st::SpaceTime, par::Collection, f::Real,
    diffop::AbstractMatrix; max_iter::Integer
) # -> Vector
    # constants
    q = @. aS - par.A + f # surface invariant flux
    diffop = get_diffop(st) # diffusion operator
    # open water
    water = vars.E .>= 0 # water set
    dw = par.cw + st.dt * par.B
    s = fill(st.dt / dw, st.nx) # diffusion coefficient
    r = (@. (vars.E + st.dt * (q + par.Fb)) / dw) # right-hand side
    # freezing ice
    h = get(vars, :h, @. -vars.E / par.Lf * (vars.E<0))
    df = @. par.k / h + par.B
    sf = @. 1 / df
    rf = @. q / df
    # interate to steady melting set
    T = get(vars, :T, @. vars.E / par.cw) # for initial guess of the melting set
    newmelt = @. T>=0 && !water # initial guess
    iter = 0 # number of iterations
    while true
        melting = copy(newmelt)
        # construct sets
        sit = copy(s)
        rit = copy(r)
        sit[melting] .= 0
        rit[melting] .= 0
        freezing = @. !(water || melting)
        sit[freezing] .= sf[freezing]
        rit[freezing] .= rf[freezing]
        # solve for T
        T .= (LA.I - par.D * SA.spdiagm(sit) * diffop) \ (0 .+ rit)
        # update melting set
        @. newmelt =  !water && q - par.B*T + par.D * $(diffop*T) >= 0
        iter += 1
        newmelt == melting && break # break if melting set is steady
        if iter > max_iter
            @warn(
                "Maximum number of iterations reached when solving for melting set. Result may be inaccurate.",
                iter, max_iter, (1:st.nx)[newmelt .!== melting]
            )
            break
        end # if >
    end # while true
    return T
end # function solveT

Infrastructure.initialise(
    model::ClassicModel, st::SpaceTime, forcing::Forcing, par::Collection, init::Collection{Vec};
    solver::AbstractSolver, lastonly::Bool, _...
) = create_storages(model, Set{Symbol}((:E, :T, :h)), st, forcing, par, init; solver, lastonly)
    # -> Tuple{Collection{Vec},Solutions{ClassicModel,F,V},Solutions{ClassicModel,F,V}}

function step_temperature!(
    ::GhostLayerSolver, vars::Collection{<:Vector}, aS::Vector, st::SpaceTime,
    par::Collection, f::Real, ::Integer, stat::NamedTuple, i::Integer
) # -> Tuple{Vector,Vector}
    C = @. aS + stat.cg_tau*vars.Tg - par.A + f
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
    return (vars.T, vars.E)
end # function specialised_step!

function step_temperature!(
    ::NonlinearSolver, vars::Collection{Vector{FT}}, aS::Vector{FT},
    st::SpaceTime, par::Collection{FT}, f::FT, _...
)::NTuple{2,Vector{FT}} where FT <: AbstractFloat
    vars.T0 = solveT0(
        get(vars, :T0, zeros(FT, st.nx)),
        vars.E, get(vars, :h, @. -vars.E / par.Lf * (vars.E<0)),
        aS, get_diffop(st), f, par
    )
    vars.T = surface_temperature(vars.E, vars.T0, par)
    C = @. aS + par.D * $(get_diffop(st)vars.T) - par.A + f
    @. vars.E += st.dt * (C - par.B*vars.T + par.Fb)
    return (vars.T, vars.E)
end # function specialised_step!

function step_temperature!(
    ::ActiveSetSolver, vars::Collection{<:Vector}, aS::Vector, st::SpaceTime,
    par::Collection, f::Real, max_iter::Integer=1000, _...
)
    diffop = get_diffop(st)
    vars.T = solveT(vars, aS, st, par, f, diffop; max_iter)
    C = @. aS + par.D * $(diffop*vars.T) - par.A + f
    @. vars.E += st.dt * (C - par.B*vars.T + par.Fb)
    return (vars.T, vars.E)
end # function specialised_step!

function Infrastructure.step!(
    ::ClassicModel, t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime, par::Collection;
    solver::AbstractSolver, _...
)::Collection{Vec}
    # get static variables
    stat = get_statics(st, par)
    # get time index
    i = round(Int, mod1((t + st.dt/2) * st.nt, st.nt))
    # forcing
    alpha = @. stat.aw * (vars.E>=0) + par.ai * (vars.E<0) # WE15 Eq. (4)
    # step T and E
    step_temperature!(solver, vars, alpha.*stat.S[:,i], st, par, f, 1000, stat, i)
    # Infer ice thickness
    vars.h = @. -vars.E / par.Lf * (vars.E<0)
    return vars
end # function Infrastructure.step!

end # module ClassicEBM
