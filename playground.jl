using EnergyBalanceModel.Infrastructure, EnergyBalanceModel.MIZEBM
using Infiltrator

# st = SpaceTime(100, 2000, 30)
st = SpaceTime(100, 10, 1)
forcing = Forcing(0.0)
params = get_defaultpar(miz_paramset)
init = Variables(
    :Ei => zeros(st.nx),
    :Ew => zeros(st.nx),
    :h => zeros(st.nx),
    :D => zeros(st.nx),
    :phi => zeros(st.nx)
)

@profview sols = integrate(st, forcing, params, init)
sols = integrate(st, forcing, params, init)
# if !(T0 isa Vector{Float64}); Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__); end

# using EnergyBalanceModel.Utilities

# to = collect(1.0:10.0)
# from = -1.0
# @code_warntype condcopy!(to, from, >(3.1))

(ice_temp(T0::VT, par::Parameters)::VT) where {VT<:AbstractVector{<:Number}} = min.(T0, par.Tm)
function T0eq(
    T0::VT,
    args::@NamedTuple{
        x::VarVec, t::Float64, h::VarVec, Tw::VarVec, phi::VarVec, f::Float64, st::SpaceTime, par::Parameters
    }
)::VT where {VT<:AbstractVector{<:Number}} # T0eq(
    conduction = @. args.par.k * (args.par.Tm - T0) / args.h
    solarin = coalbedo(args.x, true, args.par) .* solar(args.x, args.t, args.par)
    olr = @. -args.par.A - args.par.B * (T0 - args.par.Tm)
    dlap = diffusion(Tbar(ice_temp(T0, args.par), args.Tw, args.phi), args.st, args.par)
    forcing = args.f
    return @. conduction + solarin + olr + dlap + forcing
end # function T0eq


function foo(n::T)::T where {T<:Number}
    return n + 1
end

bar =
    let a = 2
        function bari(n::T)::T where {T<:Number}
            return foo(n) + a
        end

        function bar(n::T)::Nothing where {T<:Number}
            @code_warntype bari(n)
        end
    end

bar(3)
