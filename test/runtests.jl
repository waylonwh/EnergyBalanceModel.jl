using EnergyBalanceModel, Test

import EnergyBalanceModel.Infrastructure:AbstractModel

@testset "Test for annual hemispheric means" begin
    st = SpaceTime(180, 2000, 20)
    forcing = Forcing(0.0)

    mizpar = default_parameters(MIZ)
    T = fill(17.0, st.nx)
    mizinit = Collection{Vec}(
        :Ei => zeros(st.nx),
        :Ew => T .* mizpar.cw,
        :h => zeros(st.nx),
        :D => zeros(st.nx),
        :Tg => T,
    ) # Collection

    clapar = default_parameters(Classic)
    clainit = Collection{Vec}(
        :E => clapar.cw * T,
        :Tg => T
    )

    mizsols = integrate(MIZ, st, forcing, mizpar, mizinit; updatefreq=Inf)
    clasols = integrate(Classic, st, forcing, clapar, clainit; updatefreq=Inf)

    (lastyear_hemi_mean(sols::Solutions{<:AbstractModel,F,C}, var::Symbol)::Float64) where {F, C} =
        hemispheric_mean(getproperty(sols.annual.avg, var)[sols.spacetime.dur], sols.spacetime.x)

    @test lastyear_hemi_mean(mizsols, :T) - lastyear_hemi_mean(clasols, :T) < 0.1
    @test lastyear_hemi_mean(mizsols, :E) - lastyear_hemi_mean(clasols, :E) < 1.0
end # @testset begin
