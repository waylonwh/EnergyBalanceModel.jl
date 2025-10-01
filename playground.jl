using EnergyBalanceModel.Infrastructure, EnergyBalanceModel.MIZEBM

# st = SpaceTime(100, 2000, 30)
st = SpaceTime(100, 2000, 1)
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
