module ClassicEBM # EnergyBalanceModel.

using ..Infrastructure, ..Utilities

import LinearAlgebra as LA, SparseArrays as SA

@persistent(
    cg_tau = 0.0,
    dt_tau = 0.0,
    dc = 0.0,
    kappa = SA.spzeros(Float64, 0, 0),
    S = zeros(Float64, 0, 0),
    M = 0.0,
    aw = Float64[],
    kLf = 0.0,
    id::UInt = UInt(0),
    @inline function get_statics(st::SpaceTime, par::Par)
        T = eltype(st.x)
        if id != hash((st, par)) # recompute only if st or par changed
            # Difinitions for implicit scheme for Tg
            cg_tau = par.cg / par.tau
            dt_tau = st.dt / par.tau
            dc = dt_tau * cg_tau
            kappa = (one(T)+dt_tau) * LA.I(st.nx) - st.dt * par.D * get_diffop(st) / par.cg
            # Seasonal forcing [WE15 Eq. (3)]
            S = repeat(par.S0 .- par.S2 .* st.x.^2, 1, st.nt) -
                repeat(par.S1 * cos.(T(2pi) * st.t'), st.nx, 1) .* repeat(st.x, 1, st.nt)
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

Infrastructure.initialise(
    model::ClassicModel, st::SpaceTime{F,T}, forcing::Forcing{T,C}, par::Par{T}, init::Collection{Vector{T}};
    lastonly::Bool=true
) where {F,C,T<:AbstractFloat} =
    create_storages(model, Set{Symbol}((:E, :T, :h)), st, forcing, par, init; lastonly)
    # -> Tuple{Collection{Vector},Solutions{ClassicModel,F,V},Solutions{ClassicModel,F,V}}

function Infrastructure.step!(
    ::ClassicModel, t::T, f::T, vars::Collection{Vector{T}}, st::SpaceTime{F,T}, par::Par{T}
)::Collection{Vector{T}} where {F,T<:AbstractFloat}
    # get static variables
    stat = get_statics(st, par)
    # get time index
    i = round(Int, mod1((t + st.dt/2) * st.nt, st.nt))
    # forcing
    alpha = @. stat.aw * (vars.E>=0) + par.ai * (vars.E<0) # WE15 Eq. (4)
    C = @. alpha*stat.S[:,i] + stat.cg_tau*vars.Tg - par.A + f
    # surface temperature
    T0 = @. C / (stat.M - stat.kLf/vars.E) # WE15 Eq. (A3)
    vars.T = @. vars.E/par.cw * (vars.E>=0) + T0 * (vars.E<0)*(T0<0) # WE15 Eq. (9)
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
        ) # () # vars.Tg # WE15 Eq. (A1)
    # Infer ice thickness
    vars.h = @. -vars.E / par.Lf * (vars.E<0)
    return vars
end # function Infrastructure.step!

end # module ClassicEBM
