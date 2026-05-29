using EnergyBalanceModel, Test

st = SpaceTime{sin}(180, 2000, 30)
forcing = Forcing(0.0)

mizpar = default_parameters(MIZModel())
T = fill(17.0, st.nx)
mizinit = Collection{Vec}(
    :Ei => zeros(st.nx),
    :Ew => mizpar.cw * T,
    :h => zeros(st.nx),
    :D => zeros(st.nx),
    :Tg => T,
) # Collection

clapar = default_parameters(ClassicModel())
clainit = Collection{Vec}(
    :E => clapar.cw * T,
    :Tg => T
)

wimpar = default_parameters(WIModel())

lastyear_hemi_mean(sols::Solutions, var::Symbol)::Float64 =
    hemispheric_mean(getproperty(sols.annual.avg, var)[sols.spacetime.dur], sols.spacetime.x)

@testset "Code can run" begin
    global mizsols = integrate(MIZModel(), st, forcing, mizpar, mizinit; updatefreq=Inf)
    global clasols = integrate(ClassicModel(), st, forcing, clapar, clainit; updatefreq=Inf)
    global wimsols = integrate(
        WIModel(), st, forcing, wimpar, mizinit;
        spectrum=bretschneider(3.0, 9.5), updatefreq=Inf
    )
    @test mizsols isa Solutions{MIZModel,sin,false}
    @test clasols isa Solutions{ClassicModel,sin,false}
    @test wimsols isa Solutions{WIModel,sin,false}
end # @testset begin

@testset "Test for annual hemispheric means" begin
    @test lastyear_hemi_mean(mizsols, :T) - lastyear_hemi_mean(clasols, :T) < 1.0
    @test lastyear_hemi_mean(mizsols, :E) - lastyear_hemi_mean(clasols, :E) < 10.0
end # @testset begin

@testset "WIM" begin
    spectrum = bretschneider(3.0, 9.5)
    wimsols = integrate(WIModel(), st, forcing, wimpar, mizinit; spectrum, updatefreq=Inf)
    @test wimsols isa Solutions{WIModel,sin,false}
    @test EnergyBalanceModel.WIMEBM._wavenumber_ice_cache_ref[][2].key == hash(
        (Float64, spectrum.freq, wimpar.Gamma, 1e-10, wimpar)
    )
    @test wimsols.spectrum_ref[].density == spectrum.density
    @test hasproperty(wimsols.raw, :Ewave) && hasproperty(wimsols.raw, :lambda)
    # @test lastyear_hemi_mean(wimsols, :T) - lastyear_hemi_mean(mizsols, :T) < 1.0
end
