module ClassicEBM # EnergyBalanceModel.

using ..Utilities, ..Infrastructure

import LinearAlgebra as LA, SparseArrays as SA, NonlinearSolve as NlinSol

struct ClassicNonlinearModel <: AbstractModel end

const StatNT = @NamedTuple{
    cg_tau::Float64, dt_tau::Float64, dc::Float64, kappa::Matrix{Float64},
    S::Matrix{Float64}, M::Float64, aw::Vec, kLf::Float64
}

@persistent(
    cg_tau::Float64, dt_tau::Float64, dc::Float64, kappa::Matrix{Float64},
    S::Matrix{Float64}, M::Float64, aw::Vec, kLf::Float64,
    id::UInt = UInt(0),

    @inline function get_statics(st::SpaceTime{F}, par::Collection{Float64})::StatNT where F
        if id != hash((st, par)) # recompute only if st or par changed
            # Difinitions for implicit scheme for Tg
            cg_tau = par.cg / par.tau
            dt_tau = st.dt / par.tau
            dc = dt_tau * cg_tau
            kappa = (1+dt_tau) * LA.I(st.nx) - st.dt * par.D * get_diffop(st) / par.cg
            # Seasonal forcing [WE15 Eq. (3)]
            S = repeat(par.S0 .- par.S2 * st.x.^2, 1, st.nt) -
                repeat(par.S1 * cos.(2.0pi*st.t'), st.nx, 1) .* repeat(st.x, 1, st.nt)
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

function Infrastructure.initialise(
    ::ClassicModel, st::SpaceTime{F}, forcing::Forcing{C}, par::Collection{Float64}, init::Collection{Vec};
    lastonly::Bool=true
)::Tuple{Collection{Vec},Solutions{ClassicModel,F,C},Solutions{ClassicModel,F,C}} where {F, C}
    vars = deepcopy(init)
    solvars = Set{Symbol}((:E, :T, :h))
    sols = Solutions{ClassicModel}(st, forcing, par, init, solvars, lastonly)
    annusol = Solutions{ClassicModel}(st, forcing, par, init, solvars, true) # for calculating annual means
    return (vars, sols, annusol)
end # function initialise

function Infrastructure.step!(
    ::ClassicModel,
    t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime{F}, par::Collection{Float64};
)::Collection{Vec} where F
    # get static variables
    stat = get_statics(st, par)
    # get time index
    i = round(Int, mod1((t + st.dt/2.0) * st.nt, st.nt))
    # forcing
    alpha = @. stat.aw * (vars.E>=0.0) + par.ai * (vars.E<0.0) # WE15 Eq. (4)
    C = @. alpha*stat.S[:,i] + stat.cg_tau*vars.Tg - par.A + f
    # surface temperature
    T0 = @. C / (stat.M - stat.kLf/vars.E) # WE15 Eq. (A3)
    vars.T = @. vars.E/par.cw * (vars.E>=0) + T0 * (vars.E<0.0)*(T0<0.0) # WE15 Eq. (9)
    # Forward Euler for E
    @. vars.E += st.dt * (C - stat.M*vars.T + par.Fb) # WE15 Eq. (A2)
    # Implicit Euler for Tg
    vars.Tg =
        (stat.kappa - SA.spdiagm(stat.dc ./ (stat.M .- stat.kLf./vars.E) .* (T0.<0.0).*(vars.E.<0.0))) \
        (
            vars.Tg +
            (
                stat.dt_tau * (vars.E/par.cw.*(vars.E.>=0) +
                (par.ai*view(stat.S, :, i+1) .- par.A .+ f) ./ (stat.M .- stat.kLf./vars.E) .* (T0.<0.0).*(vars.E.<0.0))
            )
        ) # () # vars.Tg # WE15 Eq. (A1)
    # Infer ice thickness
    vars.h = @. -vars.E / par.Lf * (vars.E<0.0)
    return vars
end # function Infrastructure.step!

@inline function temperature(E::Vec, T0::Vector{N}, st::SpaceTime{F}, par::Collection{Float64})::Vector{N} where {N<:Number, F}
    ice = (E.<0.0)
    T = fill(zero(N), st.nx) # melting ice
    @. T[ice && T0<0.0] = T0[ice && T0<0.0] # freezing ice
    @. T[~ice] = 0.0 + E[~ice]/par.cw # open water
    return T
end # function temperature

@inline (
    T0eq(
        T0::Vector{T},
        args::@NamedTuple{
            E::Vec, h::Vec, f::Float64, i::Int, alpha::Vec, stat::StatNT, st::SpaceTime{F},
            par::Collection{Float64}
        }
    )::Vector{T}
) where {T<:Number, F} =
    @. args.par.k * (0.0-T0) / args.h  - ( # conductive flux through ice
        -args.alpha * args.stat.S[:,args.i] + # solar
        args.par.A + args.par.B * (T0-0.0) -  # OLR
        args.par.D * $(get_diffop(args.st) * temperature(args.E, T0, args.st, args.par)) - # diffusion
        args.f # forcing
    ) # (

@persistent(
    T0::Vec = zeros(Float64, 100),

    function solveT0(
        E::Vec, h::Vec, f::Float64, i::Int, alpha::Vec, stat::StatNT, st::SpaceTime{F},
        par::Collection{Float64}
    )::Vec where F
        if length(T0) != st.nx
            T0 = zeros(Float64, st.nx)
        end
        hpos = copy(h)
        condset!(hpos, 0.1, <=(0.0)) # prevent division by zero
        prob = NlinSol.NonlinearProblem(T0eq, T0, (; E, h=hpos, f, i, alpha, stat, st, par))
        sol = NlinSol.solve(prob)
        if sol.retcode != NlinSol.ReturnCode.Success
            # @warn "Nonlinear solve for T0 did not converge at $i-th step; maximum residual: $(maximum(abs.(sol.resid)))"
        end
        T0 = sol.u
        return T0 
    end # function solveT0
) # @persistent

function Infrastructure.initialise(
    ::ClassicNonlinearModel, st::SpaceTime{F}, forcing::Forcing{C}, par::Collection{Float64}, init::Collection{Vec};
    lastonly::Bool=true
)::Tuple{Collection{Vec},Solutions{ClassicNonlinearModel,F,C},Solutions{ClassicNonlinearModel,F,C}} where {F, C}
    vars = deepcopy(init)
    solvars = Set{Symbol}((:E, :T, :h))
    sols = Solutions{ClassicNonlinearModel}(st, forcing, par, init, solvars, lastonly)
    annusol = Solutions{ClassicNonlinearModel}(st, forcing, par, init, solvars, true) # for calculating annual means
    return (vars, sols, annusol)
end # function initialise

function Infrastructure.step!(
    ::ClassicNonlinearModel,
    t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime{F}, par::Collection{Float64};
)::Collection{Vec} where F
    # get static variables
    stat = get_statics(st, par)
    # get time index
    i = round(Int, mod1((t + st.dt/2.0) * st.nt, st.nt))
    # forcing
    alpha = @. stat.aw * (vars.E>=0.0) + par.ai * (vars.E<0.0) # WE15 Eq. (4)
    # Infer ice thickness
    vars.h = @. -vars.E / par.Lf * (vars.E<0.0)
    # surface temperature
    T0 = solveT0(vars.E, vars.h, f, i, alpha, stat, st, par)
    vars.T = @. vars.E/par.cw * (vars.E>=0) + T0 * (vars.E<0.0)*(T0<0.0) # WE15 Eq. (9)
    C = @. alpha * stat.S[:, i] - par.A - par.B * vars.T + par.D * $(get_diffop(st) * vars.T) + par.Fb + f
    # Forward Euler for E
    @. vars.E += st.dt * C # WE15 Eq. (A2)
    return vars
end # function Infrastructure.step!

precompile(
    Infrastructure.step!,
    (ClassicModel, Float64, Float64, Collection{Vec}, SpaceTime{identity}, Collection{Float64})
)

end # module ClassicEBM
