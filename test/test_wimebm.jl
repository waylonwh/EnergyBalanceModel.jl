# Tests for `EnergyBalanceModel.WIMEBM`, the wave-ice interaction extension of the MIZ model.

@testset "WIMEBM" begin

    par = wi_parameters()
    S = TEST_SPECTRUM
    # a coarse thickness table is enough for the unit tests; `initialise` builds a fine one
    cache = WI.cache_wavenumber!(S.freq, par, 1//100, 3)

    @testset "dispersion_relation" begin
        # with no ice and no damping the relation collapses to the deep water limit k = ω²/g
        for omega in (0.5, 1.0, 2.0)
            @test WI.dispersion_relation(ComplexF64(omega^2 / par.g), omega, 0.0, 0.0, par) ≈ 0.0 atol = 1e-12
        end # for omega
        # a non-root leaves a residual
        @test abs(WI.dispersion_relation(ComplexF64(1.0), 1.0, 0.0, 0.0, par)) > 1e-6
        # viscous damping makes the relation complex
        @test !isreal(WI.dispersion_relation(ComplexF64(0.1), 1.0, par.Gamma, 1.0, par))
    end # @testset "dispersion_relation"

    @testset "wavenumber_ice" begin
        # open water: the analytical gravity wave solution
        @test WI.wavenumber_ice(1.0, 0.0, par, 0.0) ≈ 1 / par.g
        @test WI.wavenumber_ice(2.0, 0.0, par, 0.0) ≈ 4 / par.g
        # the solved wavenumber is a root of the dispersion relation for ice too
        for h in (0.1, 1.0, 2.5)
            k = WI.wavenumber_ice(1.0, h, par, 0.0)
            @test abs(WI.dispersion_relation(k, 1.0, 0.0, h, par)) < 1e-8
            @test real(k) > 0
        end # for h
        # thicker ice lengthens the wave, so the wavenumber drops
        @test real(WI.wavenumber_ice(1.0, 3.0, par, 0.0)) < real(WI.wavenumber_ice(1.0, 0.5, par, 0.0))
        # damping introduces the attenuation rate as the imaginary part
        @test imag(WI.wavenumber_ice(1.0, 1.0, par, par.Gamma)) > 0
        @test imag(WI.wavenumber_ice(1.0, 1.0, par, 0.0)) ≈ 0 atol = 1e-12

        # the vectorised method reads the tabulated values
        kvec = WI.wavenumber_ice(S.freq, 0.0, par, 0.0)
        @test length(kvec) == length(S.freq)
        @test kvec ≈ S.freq .^ 2 ./ par.g
        @test WI.wavenumber_ice(S.freq, 1.0, par, 0.0) ≈
            WI.wavenumber_ice.(S.freq, 1.0, Ref(par), 0.0) rtol = 1e-2
        # beyond the table the wavenumbers are solved directly
        @test WI.wavenumber_ice(S.freq, 3.0, par, 0.0) ≈ WI.wavenumber_ice.(S.freq, 3.0, Ref(par), 0.0)
    end # @testset "wavenumber_ice"

    @testset "wavenumber cache" begin
        @test cache isa NTuple{2,WI.WavenumberCache{Float64}}
        @test WI._wavenumber_ice_cache_ref[] === cache
        # entry 1 holds the undamped table, entry 2 the damped one
        @test cache[1].key == hash((Float64, S.freq, 0.0, 1e-10, par))
        @test cache[2].key == hash((Float64, S.freq, par.Gamma, 1e-10, par))
        @test cache[1].dh == 0.01
        @test cache[1].hmax == 3.0
        @test size(cache[1].wavenumber) == (length(S.freq), length(0:cache[1].dh:cache[1].hmax))
        @test all(iszero, imag.(cache[1].wavenumber)) # no damping in the first table
        @test any(!iszero, imag.(cache[2].wavenumber))

        # `get_cache` picks the table by the damping parameter and reuses it
        @test WI.get_cache(S.freq, 1.0, 0.0, 1e-10, par) === cache[1]
        @test WI.get_cache(S.freq, 1.0, par.Gamma, 1e-10, par) === cache[2]
        # a different spectrum invalidates the key and forces a recompute
        @test_logs (:warn, r"cache miss") match_mode = :any WI.get_cache(
            [1.0], 1.0, 0.0, 1e-10, par
        ) # @test_logs
        @test WI._wavenumber_ice_cache_ref[] !== cache
        # restore the table the rest of the tests rely on
        WI.cache_wavenumber!(S.freq, par, 1//100, 3)
        @test WI.get_cache(S.freq, 1.0, 0.0, 1e-10, par).key == cache[1].key

        @testset "interpolate_wavenumber" begin
            @test WI.interpolate_wavenumber(0.0, cache[1]) == cache[1].wavenumber[:, 1]
            @test WI.interpolate_wavenumber(cache[1].dh, cache[1]) ≈ cache[1].wavenumber[:, 2]
            # halfway between two tabulated thicknesses is the average of the two
            @test WI.interpolate_wavenumber(cache[1].dh / 2, cache[1]) ≈
                (cache[1].wavenumber[:, 1] .+ cache[1].wavenumber[:, 2]) ./ 2
            # and the interpolant tracks the directly solved values
            @test WI.interpolate_wavenumber(0.505, cache[1]) ≈
                WI.wavenumber_ice.(S.freq, 0.505, Ref(par), 0.0) rtol = 1e-3
        end # @testset "interpolate_wavenumber"
    end # @testset "wavenumber cache"

    @testset "moments" begin
        # the zeroth moment of the spectrum gives the significant wave height, Hs = 4√m₀
        @test 4sqrt(WI.moment(S, 0)) ≈ 3.0 rtol = 0.05
        @test WI.moment_elevation(S, 0) == WI.moment(S, 0)
        @test WI.moment(S, 2) > 0
        # a constant coefficient scales the moment
        @test WI.moment(S, 0; coeff=fill(2.0, length(S.freq))) ≈ 2WI.moment(S, 0)
        @test WI.moment(S, 0; coeff=zeros(length(S.freq))) == 0.0
        # the mean wave period follows from the zeroth and second moments
        @test WI.wave_period(S) ≈ 2pi * sqrt(WI.moment(S, 0) / WI.moment(S, 2))
        @test WI.wave_period(S) > 0
    end # @testset "moments"

    @testset "strain and wavelength" begin
        strain = WI.wave_strain(S, 1.0, par)
        @test strain ≈ 2sqrt(WI.moment_strain(S, 0, 1.0, par))
        @test strain > 0
        # a calmer sea strains the ice less
        calm = bretschneider(0.5, 9.5, S.freq)
        @test WI.wave_strain(calm, 1.0, par) < strain
        # the strain is a root moment of a density that scales with Hs², so it scales with Hs
        @test WI.wave_strain(bretschneider(6.0, 9.5, S.freq), 1.0, par) ≈ 2strain

        wavelength = WI.wave_length(S, 1.0, par)
        @test wavelength isa Float64
        @test wavelength > 0
        @test wavelength ≈ 2pi / real(WI.wavenumber_ice(2pi / WI.wave_period(S), 1.0, par, 0.0)) rtol = 1e-6
        # waves are longer in thicker ice
        @test WI.wave_length(S, 3.0, par) > wavelength
    end # @testset "strain and wavelength"

    @testset "attenuation" begin
        alpha1 = WI.ice_attenuation(S, 1.0, par)
        @test length(alpha1) == length(S.freq)
        @test all(>(0), alpha1)
        @test alpha1 ≈ imag.(WI.wavenumber_ice(S.freq, 1.0, par, par.Gamma))

        atted = WI.attenuate(S, 1e5, 1.0, 0.5, par)
        @test atted isa Spectrum
        @test atted.freq == S.freq
        @test atted.period == S.period
        @test all(atted.density .< S.density)
        # nothing is attenuated over zero distance, or in the absence of ice cover
        @test WI.attenuate(S, 0.0, 1.0, 0.5, par).density ≈ S.density
        @test WI.attenuate(S, 1e5, 1.0, 0.0, par).density ≈ S.density
        # the two methods agree, and attenuation compounds over distance
        @test WI.attenuate(S, 1e5, 0.5, alpha1).density ≈ atted.density
        @test all(WI.attenuate(S, 2e5, 1.0, 0.5, par).density .< atted.density)
        # attenuation over two consecutive stretches equals attenuation over their sum
        @test WI.attenuate(WI.attenuate(S, 1e5, 0.5, alpha1), 1e5, 0.5, alpha1).density ≈
            WI.attenuate(S, 2e5, 0.5, alpha1).density
    end # @testset "attenuation"

    @testset "fracture_distance" begin
        alpha1 = WI.ice_attenuation(S, 1.0, par)
        # over this distance the strain decays from well above to well below the threshold
        L = 1e7
        @test WI.wave_strain(S, 1.0, par) > par.Ec
        @test WI.wave_strain(WI.attenuate(S, L, 0.5, alpha1), 1.0, par) < par.Ec
        l = WI.fracture_distance(S, 1.0, 0.5, L, par)
        @test 0 < l < L
        # the ice breaks up to exactly the point where the strain drops to the critical value
        @test WI.wave_strain(WI.attenuate(S, l, 0.5, alpha1), 1.0, par) ≈ par.Ec rtol = 1e-6
        # a denser ice cover attenuates faster, so the waves penetrate less far
        @test WI.fracture_distance(S, 1.0, 0.9, L, par) < l
    end # @testset "fracture_distance"

    @testset "cell_mean" begin
        alpha1 = WI.ice_attenuation(S, 1.0, par)
        # the mean of a constant is that constant
        @test WI.cell_mean(Returns(1.0), S, 0.5, alpha1, 100.0) ≈ 1.0
        @test WI.cell_mean(Returns(-2.5), S, 0.5, alpha1, 100.0) ≈ -2.5
        # the average significant wave height across a cell is below the incident one
        incident = 4sqrt(WI.moment(S, 0))
        @test WI.cell_mean(sl -> 4sqrt(WI.moment(sl, 0)), S, 0.5, alpha1, 1e5) < incident
    end # @testset "cell_mean"

    @testset "floe size" begin
        # the minimum flexural failure length grows with thickness
        @test WI.flexural_min(1.0, par) ≈
            1//2 * (pi^4 * par.Y / (48par.rhow * par.g * (1 - par.nu^2)))^(1//4)
        @test WI.flexural_min(2.0, par) > WI.flexural_min(1.0, par)
        @test WI.flexural_min(0.0, par) == 0.0

        # the truncated power law expectation lies between its bounds
        @test WI.expectation(1.0, 10.0, par.gamma) > 1.0
        @test WI.expectation(1.0, 10.0, par.gamma) < 10.0
        # a degenerate range collapses to the lower bound
        @test WI.expectation(2.0, 1.0, par.gamma) == 2.0
        @test WI.expectation(2.0, 2.0, par.gamma) == 2.0

        alpha1 = WI.ice_attenuation(S, 1.0, par)
        # the mean floe size never exceeds the size the floes already have
        @test WI.mean_size(S, 1.0, 0.5, 1e4, alpha1, par, 500.0) <= 500.0
        @test WI.mean_size(S, 1.0, 0.5, 1e4, alpha1, par, 1.0) ≈ 1.0
        @test WI.mean_size(S, 1.0, 0.5, 1e4, alpha1, par, 500.0) > 0
    end # @testset "floe size"

    @testset "grid_length" begin
        st = SpaceTime{sin}(45, 200, 1)
        lengths = WI.grid_length.(Ref(st), 1:st.nx)
        @test all(>(0), lengths)
        # the gridboxes tile a quarter of the Earth's circumference, 10⁷ m
        @test sum(lengths) ≈ 1e7
        # the sine grid is uniform in latitude, so all the gridboxes are the same length
        @test all(l -> isapprox(l, 1e7 / st.nx), lengths)

        # a grid uniform in x instead stretches towards the pole
        idlengths = WI.grid_length.(Ref(SpaceTime{identity}(45, 200, 1)), 1:45)
        @test sum(idlengths) ≈ 1e7
        @test issorted(idlengths)
        @test idlengths[end] > idlengths[1]
    end # @testset "grid_length"

    @testset "updateD!" begin
        vars = Collection{Vector{Float64}}(:D => [10.0, 20.0])
        # without a fracture distance the new size simply replaces the old one
        @test WI.updateD!(30.0, 2, vars) == 30.0
        @test vars.D == [10.0, 30.0]
        # otherwise the new size is weighted by the fraction of the cell that broke up
        @test WI.updateD!(30.0, 1, vars, 0.5, 2.0) == (0.5 * 30 + 1.5 * 10) / 2
        @test vars.D[1] == 15.0
    end # @testset "updateD!"

    @testset "initialise" begin
        st = SpaceTime{sin}(20, 100, 1)
        vars, sols, annusol = IF.initialise(
            WIModel(), st, Forcing(0.0), par, miz_initconds(st);
            solver=GhostLayerSolver(), lastonly=true, spectrum=S
        ) # IF.initialise
        @test sols isa Solutions{WIModel,sin,false}
        @test keys(sols.raw) == Set((:Ei, :Ew, :D, :h, :E, :Ti, :Tw, :T, :phi, :n))
        # the spectrum is stored by value, so later changes to the argument cannot leak in
        @test sols.spectrum_ref[].density == S.density
        @test sols.spectrum_ref[] !== S
        @test annusol.spectrum_ref[].density == S.density
        @test IF.get_spectrum(sols).freq == S.freq
        # initialisation tabulates the wavenumbers for the incident frequencies
        @test WI._wavenumber_ice_cache_ref[][2].key == hash((Float64, S.freq, par.Gamma, 1e-10, par))
        @test WI._wavenumber_ice_cache_ref[][1].hmax == 10.0
    end # @testset "initialise"

    @testset "step!" begin
        st = SpaceTime{sin}(30, 200, 3)
        vars, _, _ = IF.initialise(
            WIModel(), st, Forcing(0.0), par, miz_initconds(st);
            solver=GhostLayerSolver(), lastonly=true, spectrum=S
        ) # IF.initialise
        # with no ice anywhere the wave code is skipped and the state is untouched by it
        Dbefore = copy(vars.D)
        IF.step!(WIModel(), st.t[1], 0.0, vars, st, par; solver=GhostLayerSolver(), spectrum=S)
        @test vars.D == Dbefore

        for ti in 2:(st.nt * st.dur)
            IF.step!(
                WIModel(), st.t[mod1(ti, st.nt)], 0.0, vars, st, par;
                solver=GhostLayerSolver(), spectrum=S
            ) # IF.step!
        end # for ti
        # the same thermodynamic invariants as the MIZ model hold
        @test any(>(0), vars.phi)
        @test vars.E ≈ vars.Ei .+ vars.Ew
        @test all(<=(0), vars.Ei)
        @test all(>=(0), vars.Ew)
        @test all(p -> 0 <= p <= 1, vars.phi)
        @test all(>=(0), vars.h)
        @test all(d -> 0 <= d <= par.Dmax, vars.D)
        @test all(isfinite, vars.T)
        @test all(iszero, vars.D[iszero.(vars.Ei)])
        # waves keep the floes smaller than the same run without them
        mizvars, _, _ = IF.initialise(
            MIZModel(), st, Forcing(0.0), par, miz_initconds(st);
            solver=GhostLayerSolver(), lastonly=true, spectrum=nothing
        ) # IF.initialise
        for ti in 1:(st.nt * st.dur)
            IF.step!(
                MIZModel(), st.t[mod1(ti, st.nt)], 0.0, mizvars, st, par;
                solver=GhostLayerSolver(), spectrum=nothing
            ) # IF.step!
        end # for ti
        @test maximum(vars.D) < maximum(mizvars.D)
    end # @testset "step!"

end # @testset "WIMEBM"
