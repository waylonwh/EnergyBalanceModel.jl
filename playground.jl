using EnergyBalanceModel.Infrastructure, EnergyBalanceModel.MIZEBM
using Infiltrator
using BenchmarkTools

st = SpaceTime(100, 2000, 30)
# st = SpaceTime(100, 50, 1)
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
@benchmark sols = integrate(st, forcing, params, init)
sols = integrate(st, forcing, params, init)
# if !(T0 isa Vector{Float64}); Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__); end
