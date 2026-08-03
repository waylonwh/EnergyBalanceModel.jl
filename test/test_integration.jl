# End-to-end tests of `integrate`, `solve` and the `Solutions` they produce. These cover
# what the original `test/runtests.jl` checked, on a coarser grid and with an explicit
# solver (the default `ActiveSetSolver` is not implemented yet; see test_solvers.jl).

@testset "Integration" begin

    @testset "Code can run" begin
        @test testsol(:miz_gl) isa Solutions{MIZModel,sin,false}
        @test testsol(:classic_gl) isa Solutions{ClassicModel,sin,false}
        @test testsol(:wi_gl) isa Solutions{WIModel,sin,false}
        @test testsol(:miz_nl) isa Solutions{MIZModel,sin,false}
        @test testsol(:classic_nl) isa Solutions{ClassicModel,sin,false}
        @test testsol(:miz_var) isa Solutions{MIZModel,sin,true}
    end # @testset "Code can run"

    @testset "Test for annual hemispheric means" begin
        # the marginal ice zone changes the climate, but only mildly
        @test lastyear_hemi_mean(testsol(:miz_gl), :T) - lastyear_hemi_mean(testsol(:classic_gl), :T) < 1.0
        @test lastyear_hemi_mean(testsol(:miz_gl), :E) - lastyear_hemi_mean(testsol(:classic_gl), :E) < 10.0
        # the wave-ice interaction perturbs the MIZ climate less still
        @test lastyear_hemi_mean(testsol(:wi_gl), :T) - lastyear_hemi_mean(testsol(:miz_gl), :T) < 1.0
        # and every hemispheric mean is a sensible temperate climate
        for key in (:classic_gl, :miz_gl, :wi_gl, :classic_nl, :miz_nl)
            @test 0.0 < lastyear_hemi_mean(testsol(key), :T) < 30.0
        end # for key
    end # @testset "Test for annual hemispheric means"

    @testset "WIM" begin
        sols = testsol(:wi_gl)
        @test sols isa Solutions{WIModel,sin,false}
        # the wavenumbers were tabulated for exactly this spectrum and parameter set
        @test WI._wavenumber_ice_cache_ref[][2].key ==
            hash((Float64, TEST_SPECTRUM.freq, wi_parameters().Gamma, 1e-10, wi_parameters()))
        # and the spectrum travels with the solution
        @test sols.spectrum_ref[].density == TEST_SPECTRUM.density
        @test sols.spectrum_ref[].freq == TEST_SPECTRUM.freq
        # waves break the floes up, so they stay smaller than in the MIZ model alone
        @test maximum(maximum, sols.annual.avg.D) < maximum(maximum, testsol(:miz_gl).annual.avg.D)
        # a spectrum is compulsory for the wave-ice model
        @test_throws ArgumentError integrate(
            WIModel(), SpaceTime{sin}(10, 50, 1), Forcing(0.0), wi_parameters(),
            miz_initconds(SpaceTime{sin}(10, 50, 1)); solver=GhostLayerSolver(), updatefreq=Inf
        ) # @test_throws
    end # @testset "WIM"

    @testset "Solutions contents" begin
        sols = testsol(:miz_gl)
        @test sols.spacetime === ST
        @test sols.forcing isa Forcing{false}
        @test keys(sols.raw) == Set((:Ei, :Ew, :D, :h, :E, :Ti, :Tw, :T, :phi, :n))
        @test keys(sols.initconds) == Set((:Tg, :Ei, :Ew, :h, :D))
        @test length(sols.ts) == ST.nt
        @test all(var -> length(sols.raw[var]) == ST.nt, keys(sols.raw))
        @test all(var -> length(sols.annual.avg[var]) == ST.dur, keys(sols.raw))
        @test all(snapshot -> length(snapshot) == ST.nx, sols.raw.T)
        # every stored field is finite, apart from the ice temperature where there is no ice
        @test all(snapshot -> all(isfinite, snapshot), sols.raw.T)
        @test all(snapshot -> all(isfinite, snapshot), sols.raw.E)
        @test all(snapshot -> all(isfinite, snapshot), sols.annual.avg.phi)
        # the annual average is the mean of the stored year
        @test sols.annual.avg.T[ST.dur] ≈ UT.crossmean(sols.raw.T)
        # the seasonal snapshots are the corresponding raw timesteps of the last year
        @test sols.annual.winter.T[ST.dur] == sols.raw.T[ST.winter.inx]
        @test sols.annual.summer.T[ST.dur] == sols.raw.T[ST.summer.inx]
        # ... and summer really is warmer than winter
        @test hemispheric_mean(sols.annual.summer.T[ST.dur], ST.x) >
            hemispheric_mean(sols.annual.winter.T[ST.dur], ST.x)

        @testset "show" begin
            @test occursin("Solutions{MIZModel, sin, false}", repr(sols))
            pretty = sprint(show, MIME"text/plain"(), sols)
            @test occursin("10 solution variables", pretty)
            @test occursin("45 latitudinal gridboxes", pretty)
            @test occursin("200 timesteps", pretty)
            @test occursin("Forcing(0.0)", pretty)
        end # @testset "show"
    end # @testset "Solutions contents"

    @testset "lastonly" begin
        full = testsol(:classic_full)
        @test !full.lastonly
        @test length(full.ts) == ST.nt * ST.dur
        @test full.ts[1] ≈ ST.dt / 2
        @test full.ts[end] ≈ ST.dur - ST.dt / 2
        @test length(full.raw.T) == ST.nt * ST.dur
        # the last year of the full record is exactly what the default run stores
        @test full.raw.T[(end - ST.nt + 1):end] == testsol(:classic_gl).raw.T
        @test full.annual.avg.T == testsol(:classic_gl).annual.avg.T
    end # @testset "lastonly"

    @testset "ice_area from a Solutions" begin
        sols = testsol(:miz_gl)
        for season in (:winter, :summer, :avg)
            area = ice_area(sols, season, ST.dur)
            @test area >= 0
            @test area <= ice_area(ones(ST.nx), ST.x)
        end # for season
        # the hemisphere carries most of its ice in winter
        @test ice_area(sols, :winter, ST.dur) >= ice_area(sols, :summer, ST.dur)
        @test ice_area(sols, :avg, ST.dur) ≈ ice_area(sols.annual.avg.phi[ST.dur], ST.x)
        # the ClassicModel method passes a BitVector to a method that only takes Vec
        @test_broken ice_area(testsol(:classic_gl), :winter, ST.dur) isa Float64
    end # @testset "ice_area from a Solutions"

    @testset "soldiff" begin
        miz = testsol(:miz_gl)
        classic = testsol(:classic_gl)
        diff = soldiff(miz, classic)
        @test diff isa Solutions{IF.ModelDiff{MIZModel,ClassicModel},sin,false}
        @test (miz - classic).raw.T[1] == diff.raw.T[1] # `soldiff` is the `-` operator
        @test diff.spacetime === miz.spacetime
        # only the variables the two models share survive
        @test keys(diff.raw) == intersect(keys(miz.raw), keys(classic.raw))
        @test keys(diff.raw) == Set((:E, :T, :h))
        # and the entries really are differences
        @test diff.raw.T[1] == miz.raw.T[1] .- classic.raw.T[1]
        @test diff.raw.E[end] == miz.raw.E[end] .- classic.raw.E[end]
        @test diff.annual.avg.T[ST.dur] == miz.annual.avg.T[ST.dur] .- classic.annual.avg.T[ST.dur]
        @test diff.annual.winter.T[1] == miz.annual.winter.T[1] .- classic.annual.winter.T[1]
        # the forcings subtract too
        @test diff.forcing.base == miz.forcing.base - classic.forcing.base
        # a model differenced with itself is identically zero
        selfdiff = soldiff(miz, miz)
        @test all(iszero, selfdiff.raw.T[1])
        @test all(iszero, selfdiff.annual.avg.E[ST.dur])

        # solutions on different grids cannot be compared
        otherst = SpaceTime{sin}(20, 200, 3)
        other = integrate(
            ClassicModel(), otherst, Forcing(0.0), classic_parameters(), classic_initconds(otherst);
            solver=GhostLayerSolver(), updatefreq=Inf
        ) # integrate
        @test_throws ArgumentError soldiff(miz, other)
        # nor can solutions under a varying forcing, whose difference is not a `Forcing`
        @test_throws MethodError soldiff(testsol(:miz_var), testsol(:classic_var))
    end # @testset "soldiff"

    @testset "solve" begin
        prob = EBMProblem(
            ClassicModel();
            st=(; nx=20, nt=100, dur=2), forcing=0.0, initconds=fill(17.0, 20)
        ) # EBMProblem
        sols = solve(prob, GhostLayerSolver(); updatefreq=Inf)
        @test sols isa Solutions{ClassicModel,sin,false}
        @test sols.spacetime === prob.st
        @test sols.parameters === prob.parameters
        @test sols.initconds === prob.initconds
        @test sols.solver isa GhostLayerSolver
        @test sols.lastonly
        @test !solve(prob, GhostLayerSolver(); lastonly=false, updatefreq=Inf).lastonly
        # `solve` only unpacks the problem and calls `integrate`
        direct = integrate(
            prob.model, prob.st, prob.forcing, prob.parameters, prob.initconds;
            solver=GhostLayerSolver(), updatefreq=Inf
        ) # integrate
        @test sols.raw.T[end] == direct.raw.T[end]
        # a `WIModel` problem carries its spectrum through to `integrate`
        wiprob = EBMProblem(
            WIModel(); st=(; nx=10, nt=50, dur=1), spectrum=TEST_SPECTRUM
        ) # EBMProblem
        wisols = solve(wiprob, GhostLayerSolver(); updatefreq=Inf)
        @test wisols isa Solutions{WIModel,sin,false}
        @test wisols.spectrum_ref[].density == TEST_SPECTRUM.density
    end # @testset "solve"

    @testset "progress reporting" begin
        st = SpaceTime{sin}(10, 50, 1)
        out = capture_stdout() do
            integrate(
                ClassicModel(), st, Forcing(0.0), classic_parameters(), classic_initconds(st);
                solver=GhostLayerSolver(), updatefreq=0.0
            ) # integrate
        end # do
        @test occursin("Integrating ClassicModel", out)
        @test occursin("50/50", out)
        @test occursin("Done", out)
        @test occursin("t = ", out) # the custom info feed
        # an infinite update frequency silences the progress bar
        quiet = capture_stdout() do
            integrate(
                ClassicModel(), st, Forcing(0.0), classic_parameters(), classic_initconds(st);
                solver=GhostLayerSolver(), updatefreq=Inf
            ) # integrate
        end # do
        @test isempty(quiet)
    end # @testset "progress reporting"

    @testset "climate response to forcing" begin
        # warming the climate must warm it: compare the last year of two constant forcings
        st = SpaceTime{sin}(30, 200, 3)
        cold = integrate(
            ClassicModel(), st, Forcing(0.0), classic_parameters(), classic_initconds(st);
            solver=GhostLayerSolver(), updatefreq=Inf
        ) # integrate
        warm = integrate(
            ClassicModel(), st, Forcing(10.0), classic_parameters(), classic_initconds(st);
            solver=GhostLayerSolver(), updatefreq=Inf
        ) # integrate
        @test lastyear_hemi_mean(warm, :T) > lastyear_hemi_mean(cold, :T)
        @test ice_area((warm.annual.winter.E[3] .< 0) .* 1.0, st.x) <
            ice_area((cold.annual.winter.E[3] .< 0) .* 1.0, st.x)
    end # @testset "climate response to forcing"

    @testset "run_example" begin
        # too slow to run here (50 years on a 180-point grid), so only check the interface
        @test hasmethod(run_example, Tuple{ClassicModel})
        @test hasmethod(run_example, Tuple{MIZModel})
        @test hasmethod(run_example, Tuple{WIModel})
        @test hasmethod(run_example, Tuple{})
    end # @testset "run_example"

end # @testset "Integration"
