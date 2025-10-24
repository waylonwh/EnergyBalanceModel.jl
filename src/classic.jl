module Classic # EnergyBalanceModel.

using ..Utilities, ..Infrastructure

import LinearAlgebra as LA, SparseArrays as SA

@persistent(
    cg_tau::Float64, dt_tau::Float64, dc::Float64, kappa::Matrix{Float64},
    S::Matrix{Float64}, M::Float64, aw::Vec, kLf::Float64,
    id::UInt64 = UInt64(0),

    @inline function get_statics(st::SpaceTime{F}, par::Collection{Float64})::@NamedTuple{
        cg_tau::Float64, dt_tau::Float64, dc::Float64, kappa::Matrix{Float64},
        S::Matrix{Float64}, M::Float64, aw::Vec, kLf::Float64
    } where F
        if id != hash((st, par)) # recompute only if st or par changed
            # Difinitions for implicit scheme for Tg
            cg_tau = par.cg / par.tau
            dt_tau = st.dt / par.tau
            dc = dt_tau * cg_tau
            kappa = (1+dt_tau) * LA.I(st.nx) - st.dt * par.D * get_diffop(st.nx) / par.cg
            # Seasonal forcing [WE15 Eq. (3)]
            S = repeat(par.S0 .- par.S2 * st.x.^2, 1, st.nt) -
                repeat(par.S1 * cos.(2.0*pi*st.t'), st.nx, 1) .* repeat(st.x, 1, st.nt)
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

function Infrastructure.step!(
    ::Val{:Classic},
    t::Float64, f::Float64, vars::Collection{Vec}, st::SpaceTime{F}, par::Collection{Float64};
    debug::Union{Expr,Nothing}=nothing
)::Collection{Vec} where F
    # get static variables
    stat = get_statics(st, par)
    # get time index
    i = round(Int, mod1((t + st.dt/2.0) * st.nt, st.nt))
    # forcing
    alpha = @. stat.aw * (vars.E>0.0) + par.ai * (vars.E<0.0) # WE15 Eq. (4)
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
    # debug
    if !isnothing(debug)
        vars.debug = eval(debug) # !
    end # if !=
    return vars
end # function Infrastructure.step!

precompile(
    Infrastructure.step!,
    (Val{:Classic}, Float64, Float64, Collection{Vec}, SpaceTime{identity}, Collection{Float64})
)

end # module Classic
