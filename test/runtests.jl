using EnergyBalanceModel, Test

st = SpaceTime{sin}(180, 2000, 30)
forcing = Forcing(0.0)

FT = eltype(st.x)
mizpar = default_parameters(MIZModel(), FT)
T = fill(FT(17.0), st.nx)
mizinit = Collection{Vector{FT}}(
    :Ei => zeros(FT, st.nx),
    :Ew => mizpar.cw * T,
    :h => zeros(FT, st.nx),
    :D => zeros(FT, st.nx),
    :Tg => T,
) # Collection

clapar = default_parameters(ClassicModel(), FT)
clainit = Collection{Vector{FT}}(
    :E => clapar.cw * T,
    :Tg => T
)

wimpar = default_parameters(WIModel(), FT)

lastyear_hemi_mean(sols::Solutions, var::Symbol) =
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

@testset "Float32 compatibility" begin
    st32 = SpaceTime{sin,Float32}(90, 500, 5)
    forcing32 = Forcing(Float32(0))
    par32 = default_parameters(MIZModel(), Float32)
    T32 = fill(Float32(5), st32.nx)
    init32 = Collection{Vector{Float32}}(
        :Ei => zeros(Float32, st32.nx),
        :Ew => par32.cw * T32,
        :h => zeros(Float32, st32.nx),
        :D => zeros(Float32, st32.nx),
        :Tg => T32,
    )
    sol32 = integrate(MIZModel(), st32, forcing32, par32, init32; updatefreq=Inf)
    @test eltype(sol32.spacetime.x) === Float32
    @test eltype(sol32.raw.T[end]) === Float32
end

# TODO add tests for WIM
