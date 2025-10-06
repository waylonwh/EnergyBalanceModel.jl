using EnergyBalanceModel.Infrastructure, EnergyBalanceModel.MIZEBM
# using Infiltrator
# Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__)

# st = SpaceTime(100, 2000, 1)
st = SpaceTime(100, 2000, 30)
# st = SpaceTime(100, 2000, 200)
forcing = Forcing(0.0)
params = get_defaultpar(miz_paramset)
init = Collection{Vec}(
    :Ei => zeros(st.nx),
    :Ew => zeros(st.nx),
    :h => zeros(st.nx),
    :D => zeros(st.nx),
    :phi => zeros(st.nx)
)

sols = integrate(st, forcing, params, init);

@profview sols = integrate(st, forcing, params, init)

# using BenchmarkTools
# @benchmark sols = integrate(st, forcing, params, init)

# using GLMakie
# fig = contourf(reduce(hcat, sols.raw.Ei)')
