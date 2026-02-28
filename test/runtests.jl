using EnergyBalanceModel, Test

import EnergyBalanceModel.Infrastructure: ClassicModel, MIZModel

st = SpaceTime{sin}(180, 2000, 30)
forcing = Forcing(0.0)

mizpar = default_parameters(miz)
T = fill(17.0, st.nx)
mizinit = Collection{Vec}(
    :Ei => zeros(st.nx),
    :Ew => mizpar.cw * T,
    :h => zeros(st.nx),
    :D => zeros(st.nx),
    :Tg => T,
) # Collection

clapar = default_parameters(classic)
clainit = Collection{Vec}(
    :E => clapar.cw * T,
    :Tg => T
)

lastyear_hemi_mean(sols::Solutions, var::Symbol)::Float64 =
    hemispheric_mean(getproperty(sols.annual.avg, var)[sols.spacetime.dur], sols.spacetime.x)

@testset "Code can run" begin
    global mizsols = integrate(miz, st, forcing, mizpar, mizinit; updatefreq=Inf)
    global clasols = integrate(classic, st, forcing, clapar, clainit; updatefreq=Inf)
    @test mizsols isa Solutions{MIZModel,sin,true}
    @test clasols isa Solutions{ClassicModel,sin,true}
end # @testset begin

@testset "Test for annual hemispheric means" begin
    @test lastyear_hemi_mean(mizsols, :T) - lastyear_hemi_mean(clasols, :T) < 1.0
    @test lastyear_hemi_mean(mizsols, :E) - lastyear_hemi_mean(clasols, :E) < 10.0
end # @testset begin

@testset "Check the solutions include T0" begin
    @test hasproperty(mizsols.raw, :T0)
    @test hasproperty(clasols.raw, :T0)
end # @testset begin
