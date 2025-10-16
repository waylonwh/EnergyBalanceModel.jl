using EnergyBalanceModel, Test, JLD2

#=
st = SpaceTime{sin}(180, 2000, 1)
forcing = Forcing(0.0)
par = default_parameters(:MIZ)
init = Collection{Vec}(
           :Ei => zeros(st.nx),
           :Ew => zeros(st.nx),
           :h => zeros(st.nx),
           :D => zeros(st.nx),
           :phi => zeros(st.nx)
       )
sols = integrate(:MIZ, st, forcing, par, init)

using JLD2
jldsave("solution_1year.jld2", sols=sols)
=#

@testset "Basic Test 1year solution..." begin

    st = SpaceTime{sin}(180, 2000, 1)
    forcing = Forcing(0.0)
    par = default_parameters(:MIZ)
    init = Collection{Vec}(
            :Ei => zeros(st.nx),
            :Ew => zeros(st.nx),
            :h => zeros(st.nx),
            :D => zeros(st.nx),
            :phi => zeros(st.nx)
        )
    sols = integrate(:MIZ, st, forcing, par, init)

    file = jldopen("solution_1year.jld2")
    sols_loaded = file["sols"]

    for key in keys(sols.raw)

        @show key
        solv = getproperty(sols.raw, key)[10]
        solv_loaded = getproperty(sols_loaded.raw, key)[10]
        EnergyBalanceModel.Utilities.condset!(solv, 0.0, isnan)
        EnergyBalanceModel.Utilities.condset!(solv_loaded, 0.0, isnan)

        @show all(isapprox.(solv, solv_loaded))
        @test all(isapprox.(solv, solv_loaded))
    end
end
