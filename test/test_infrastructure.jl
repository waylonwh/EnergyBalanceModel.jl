# Tests for the types and helpers in `EnergyBalanceModel.Infrastructure`.

@testset "Infrastructure" begin

    @testset "Exported names" begin
        exported = names(EnergyBalanceModel)
        for name in (
            :ClassicModel, :MIZModel, :WIModel,
            :ActiveSetSolver, :GhostLayerSolver, :NonlinearSolver,
            :Collection, :EBMProblem, :Forcing, :Solutions, :SpaceTime,
            :Spectrum, :bretschneider, :monochromatic,
            :default_parameters, :integrate, :soldiff, :solve,
            :hemispheric_mean, :ice_area,
            :Layout, :backend, :default_layout, :plot_avg, :plot_raw, :plot_seasonal,
            :run_example
        )
            @test name in exported
        end # for name
    end # @testset "Exported names"

    @testset "Model and solver types" begin
        @test ClassicModel() isa IF.AbstractModel
        @test MIZModel() isa IF.AbstractModel
        @test WIModel() isa IF.AbstractModel
        @test IF.ModelDiff{MIZModel,ClassicModel}() isa IF.AbstractModel
        @test ActiveSetSolver() isa IF.AbstractSolver
        @test GhostLayerSolver() isa IF.AbstractSolver
        @test NonlinearSolver() isa IF.AbstractSolver
        @test IF.DiffSolver{GhostLayerSolver,NonlinearSolver}() isa IF.AbstractSolver
        # the models and solvers are singletons
        @test MIZModel() === MIZModel()
        @test ActiveSetSolver() === ActiveSetSolver()
        @test IF.Vec === Vector{Float64}
    end # @testset "Model and solver types"

    @testset "Collection" begin
        coll = Collection{Float64}(:D => 0.6, :A => 193.0, :B => 2.1)
        @test coll isa AbstractDict{Symbol,Float64}
        @test length(coll) == 3
        @test coll.D == 0.6
        @test coll[:A] == 193.0
        @test getproperty(coll, :B) == 2.1
        @test propertynames(coll) == Set((:D, :A, :B))
        @test keys(coll) == Set((:D, :A, :B))
        @test sort(collect(values(coll))) == [0.6, 2.1, 193.0]
        @test haskey(coll, :D) && !haskey(coll, :nope)
        @test get(coll, :D, -1.0) == 0.6
        @test get(coll, :nope, -1.0) == -1.0

        # setting through both syntaxes
        coll.F = 0.0
        @test coll.F == 0.0
        coll[:G] = 1.0
        @test coll.G == 1.0
        @test length(coll) == 5

        # iteration yields Symbol => value pairs
        @test sort(collect(coll), by=first) == sort(
            [:D => 0.6, :A => 193.0, :B => 2.1, :F => 0.0, :G => 1.0], by=first
        )
        # equal contents hash equally, which `ClassicEBM.get_statics` relies on
        @test hash(Collection{Float64}(:a => 1.0)) == hash(Collection{Float64}(:a => 1.0))
        @test hash(Collection{Float64}(:a => 1.0)) != hash(Collection{Float64}(:a => 2.0))
        # an empty Collection is valid, and other element types work too
        @test isempty(Collection{Float64}())
        @test Collection{Vector{Float64}}(:x => [1.0, 2.0]).x == [1.0, 2.0]
        @test_throws KeyError Collection{Float64}(:a => 1.0).missingkey
    end # @testset "Collection"

    @testset "uniqueunion" begin
        left = Collection{Float64}(:x => 1.0, :y => 2.0)
        # when the two collections agree on every shared key nothing is recorded
        @test isempty(IF.uniqueunion(left, Collection{Float64}(:x => 1.0, :z => 3.0)))
        # a disagreement keeps the union of all keys
        merged = IF.uniqueunion(left, Collection{Float64}(:x => 9.0, :z => 3.0))
        @test keys(merged) == Set((:x, :y, :z))
        # disjoint collections trivially "agree"
        @test isempty(IF.uniqueunion(left, Collection{Float64}(:z => 3.0)))
        # the value type is widened to hold both
        @test IF.uniqueunion(left, Collection{Int}(:x => 1)) isa Collection{Real}
    end # @testset "uniqueunion"

    @testset "SpaceTime" begin
        st = SpaceTime{sin}(180, 2000, 30)
        @test st isa SpaceTime{sin}
        @test st.nx == 180
        @test st.nt == 2000
        @test st.dur == 30
        @test st.dt == 1 / 2000
        @test length(st.x) == 180
        @test issorted(st.x)
        @test 0 < st.x[1] && st.x[end] < 1 # cell centres, singularities excluded
        @test st.x ≈ sin.(st.u)
        @test length(st.t) == 2000
        @test st.t[1] ≈ st.dt / 2
        @test st.t[end] ≈ 1 - st.dt / 2
        @test length(st.T) == 2000 * 30
        @test st.T[1] ≈ st.dt / 2
        @test st.T[end] ≈ 30 - st.dt / 2
        @test st.winter == (t=0.26125, inx=round(Int, 2000 * 0.26125))
        @test st.summer == (t=0.77375, inx=round(Int, 2000 * 0.77375))

        # the default map is the identity, with a uniform grid on [0,1]
        sti = SpaceTime(10, 50, 2)
        @test sti isa SpaceTime{identity}
        @test sti.x ≈ 0.05:0.1:0.95
        @test SpaceTime{identity}(10, 50, 2).x == sti.x
        # explicit ranges and custom seasonal peaks
        stc = SpaceTime{identity}((0.0, 2.0), 4, 10, 1; winter=0.1, summer=0.6)
        @test stc.x ≈ [0.25, 0.75, 1.25, 1.75]
        @test stc.winter == (t=0.1, inx=1)
        @test stc.summer == (t=0.6, inx=6)

        @test repr(st) == "SpaceTime{sin}(180, 2000, 30)"
        pretty = sprint(show, MIME"text/plain"(), st)
        @test occursin("180 latitudinal gridboxes", pretty)
        @test occursin("2000 timesteps per year", pretty)
        @test occursin("30 years of simulation", pretty)
        @test occursin("winter at t=0.26125", pretty)
    end # @testset "SpaceTime"

    @testset "Forcing" begin
        @testset "constant" begin
            forcing = Forcing(2.5)
            @test forcing isa Forcing{false}
            @test forcing.base == 2.5
            @test forcing(0.0) == 2.5
            @test forcing(1e6) == 2.5
            @test forcing.domain == (0, 0, 0, 0, 0)
            @test repr(forcing) == "Forcing(2.5) (constant forcing)"
            @test occursin("is constant", sprint(show, MIME"text/plain"(), forcing))
        end # @testset "constant"

        @testset "variable" begin
            forcing = Forcing(0.0, 5.0, -5.0, (10, 10), (0.5, -0.5))
            @test forcing isa Forcing{true}
            @test forcing.domain == (0, 10, 20, 30, 50)
            @test forcing(0.0) == 0.0 # base
            @test forcing(9.99) == 0.0
            @test forcing(15.0) == 2.5 # warming
            @test forcing(17.57) ≈ 3.785
            @test forcing(20.0) == 5.0 # peak
            @test forcing(29.0) == 5.0
            @test forcing(40.0) == 0.0 # cooling
            @test forcing(50.0) == -5.0 # cool
            @test forcing(1e6) == -5.0
            # continuity across the four internal breakpoints
            for t in (10.0, 20.0, 30.0, 50.0)
                @test forcing(t) ≈ forcing(t - 1e-9) atol = 1e-6
            end # for t
            @test repr(forcing) == "Forcing(0.0, 5.0, -5.0)"
            pretty = sprint(show, MIME"text/plain"(), forcing)
            @test occursin("(warming)", pretty)
            @test occursin("(cooling)", pretty)

            # the warming and cooling durations must be positive whole numbers of years
            @test_throws ArgumentError Forcing(0.0, 5.0, -5.0, (10, 10), (-0.5, -0.5))
            @test_throws ArgumentError Forcing(0.0, 5.0, -5.0, (10, 10), (0.3, -0.5))
            @test_throws ArgumentError Forcing(0.0, 5.0, -5.0, (10, 10), (0.5, 0.5))
            @test_throws ArgumentError Forcing(0.0, 5.0, -5.0, (10, 10), (0.5, -0.3))
        end # @testset "variable"

        @testset "annual_mean" begin
            st = SpaceTime{sin}(180, 2000, 30)
            forcing = Forcing(0.0, 10.0, -5.0, (20, 10), (0.5, -0.5))
            @test IF.annual_mean(forcing, st, 24) ≈ 1.75
            @test IF.annual_mean(forcing, st, 1) == 0.0
            @test IF.annual_mean(Forcing(3.0), st, 7) == 3.0
        end # @testset "annual_mean"
    end # @testset "Forcing"

    @testset "Spectrum" begin
        spectrum = Spectrum([0.5, 1.0, 1.5], [0.2, 0.5, 0.1])
        @test spectrum.freq == [0.5, 1.0, 1.5]
        @test spectrum.period ≈ 2pi ./ [0.5, 1.0, 1.5]
        @test spectrum.density == [0.2, 0.5, 0.1]
        @test_throws ArgumentError Spectrum([0.5, 1.0], [0.2])
        @test repr(spectrum) == "Spectrum(3 components)"
        @test occursin("3 frequency components", sprint(show, MIME"text/plain"(), spectrum))

        @testset "bretschneider" begin
            S = bretschneider(3.0, 9.5)
            @test S isa Spectrum
            @test S.freq == collect(range(2pi/23.8, 2pi/2.5; step=7.5e-2))
            @test length(S.freq) == length(S.density) == length(S.period)
            @test all(>(0), S.density)
            # the analytical peak of the Bretschneider spectrum sits at T = Tp
            @test S.period[argmax(S.density)] ≈ 9.5 rtol = 0.05
            # the zeroth moment of the spectrum recovers Hs = 4√m₀
            @test 4sqrt(WI.moment(S, 0)) ≈ 3.0 rtol = 0.05
            # a larger significant wave height scales the density by Hs²
            @test bretschneider(6.0, 9.5).density ≈ 4S.density
            # the frequencies may be chosen freely
            @test bretschneider(3.0, 9.5, [1.0, 2.0]).freq == [1.0, 2.0]
        end # @testset "bretschneider"

        @testset "monochromatic" begin
            M = monochromatic(3.0, 9.5)
            @test M isa Spectrum
            @test M.freq[argmax(M.density)] ≈ 2pi / 9.5 atol = 1e-3
            @test 4sqrt(WI.moment(M, 0)) ≈ 3.0 rtol = 0.05
            # a narrower Gaussian gives a sharper peak but the same total energy
            @test maximum(monochromatic(3.0, 9.5; eps=1e-7).density) > maximum(M.density)
        end # @testset "monochromatic"
    end # @testset "Spectrum"

    @testset "default_parameters" begin
        classic = default_parameters(ClassicModel())
        miz = default_parameters(MIZModel())
        wim = default_parameters(WIModel())
        @test classic isa Collection{Float64}
        @test keys(classic) == Set(
            (:D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :cg, :tau)
        )
        # each model extends the parameter set of the simpler one
        @test issubset(keys(classic), keys(miz))
        @test issubset(keys(miz), keys(wim))
        @test setdiff(keys(miz), keys(classic)) ==
            Set((:Tm, :m1, :m2, :alpha, :rl, :Dmin, :Dmax, :hmin, :kappa, :rhow, :cp, :ch, :u0))
        @test setdiff(keys(wim), keys(miz)) == Set((:Y, :nu, :g, :Ec, :Gamma, :gamma))
        @test keys(wim) == keys(IF.default_parval)
        # the values come from the documented defaults
        @test classic.D == 0.6
        @test classic.A == 193.0
        @test classic.B == 2.1
        @test classic.cw == 9.8
        @test miz.Tm == 0.0
        @test miz.Dmax == 500.0
        @test wim.nu == 0.3
        @test wim.g == 9.81
        # the seconds-per-year conversions
        @test miz.m1 ≈ 1.6e-6 * 31536000
        @test miz.cp ≈ 3980.0 / 31536000
        @test wim.gamma == 2 + log2(0.9)
        # each call returns a fresh Collection, so mutation cannot leak into the defaults
        mutated = default_parameters(ClassicModel())
        mutated.D = 99.0
        @test default_parameters(ClassicModel()).D == 0.6
        @test IF.default_parval.D == 0.6
    end # @testset "default_parameters"

    @testset "initconds_fromT and check_initconds" begin
        T = fill(17.0, 5)
        classic = IF.initconds_fromT(ClassicModel(), T, 9.8)
        @test keys(classic) == Set((:Tg, :E))
        @test classic.E ≈ 9.8T
        @test classic.Tg == T

        miz = IF.initconds_fromT(MIZModel(), T, 9.8)
        @test keys(miz) == Set((:Tg, :Ei, :Ew, :h, :D))
        @test miz.Ew ≈ 9.8T
        @test miz.Ei == zeros(5)
        @test miz.h == zeros(5)
        @test miz.D == zeros(5)
        @test IF.initconds_fromT(WIModel(), T, 9.8) |> keys == keys(miz)

        @test IF.initvarset(ClassicModel()) == Set((:Tg, :E))
        @test IF.initvarset(MIZModel()) == Set((:Tg, :Ei, :Ew, :h, :D))

        @test isnothing(IF.check_initconds(ClassicModel(), classic))
        @test isnothing(IF.check_initconds(MIZModel(), miz))
        @test_throws ArgumentError IF.check_initconds(MIZModel(), classic)
        @test_throws ArgumentError IF.check_initconds(ClassicModel(), Collection{Vector{Float64}}(:E => T))
        # extra variables are tolerated, but warned about
        extra = Collection{Vector{Float64}}(:Tg => T, :E => T, :spare => T)
        @test_logs (:warn, r"only needs initial conditions") IF.check_initconds(ClassicModel(), extra)
    end # @testset "initconds_fromT and check_initconds"

    @testset "EBMProblem" begin
        @testset "defaults" begin
            prob = EBMProblem(ClassicModel())
            @test prob isa EBMProblem{Float64}
            @test prob.model isa ClassicModel
            @test prob.st isa SpaceTime{sin}
            @test (prob.st.nx, prob.st.nt, prob.st.dur) == (180, 2000, 50)
            @test prob.forcing isa Forcing{false}
            @test prob.forcing.base == 0.0
            @test keys(prob.parameters) == keys(default_parameters(ClassicModel()))
            @test keys(prob.initconds) == Set((:Tg, :E))
            @test all(prob.initconds.Tg .== 17.0)
            @test isnothing(prob.spectrum)
        end # @testset "defaults"

        @testset "keyword overrides" begin
            prob = EBMProblem(MIZModel(); st=(; nx=90, nt=500, dur=7, F=identity))
            @test prob.st isa SpaceTime{identity}
            @test (prob.st.nx, prob.st.nt, prob.st.dur) == (90, 500, 7)
            @test length(prob.initconds.Tg) == 90

            @test EBMProblem(ClassicModel(); st=ST).st === ST
            @test EBMProblem(ClassicModel(); forcing=2.0).forcing.base == 2.0
            @test EBMProblem(ClassicModel(); forcing=Forcing(3.0)).forcing.base == 3.0
            varying = EBMProblem(
                ClassicModel();
                forcing=(base=0.0, peak=5.0, cool=-5.0, holdyrs=(10, 10), rates=(0.5, -0.5))
            ).forcing
            @test varying isa Forcing{true}
            @test varying.domain == (0, 10, 20, 30, 50)

            # a NamedTuple of parameters overrides the model defaults entry by entry
            pars = EBMProblem(ClassicModel(); parameters=(; D=1.5, A=100.0)).parameters
            @test pars.D == 1.5
            @test pars.A == 100.0
            @test pars.B == default_parameters(ClassicModel()).B
            @test keys(pars) == keys(default_parameters(ClassicModel()))
            @test EBMProblem(ClassicModel(); parameters=Collection{Float64}(:cw => 1.0)).parameters.cw == 1.0

            # a temperature vector is expanded into the full initial condition set
            init = EBMProblem(MIZModel(); st=(; nx=5), initconds=fill(3.0, 5)).initconds
            @test keys(init) == Set((:Tg, :Ei, :Ew, :h, :D))
            @test init.Tg == fill(3.0, 5)
            @test init.Ew ≈ default_parameters(MIZModel()).cw * fill(3.0, 5)
            given = Collection{Vector{Float64}}(:E => zeros(5), :Tg => zeros(5))
            @test EBMProblem(ClassicModel(); st=(; nx=5), initconds=given).initconds === given
        end # @testset "keyword overrides"

        @testset "spectrum handling" begin
            @test EBMProblem(WIModel()).spectrum isa Spectrum
            @test EBMProblem(WIModel()).spectrum.density ≈ bretschneider(3.0, 9.5).density
            @test EBMProblem(WIModel(); spectrum=TEST_SPECTRUM).spectrum === TEST_SPECTRUM
            @test_logs (:warn, r"ignored for non-WIModel") EBMProblem(
                ClassicModel(); st=(; nx=5), spectrum=TEST_SPECTRUM
            )
        end # @testset "spectrum handling"

        @testset "invalid input" begin
            @test_throws ArgumentError EBMProblem(
                MIZModel(); st=(; nx=5),
                initconds=Collection{Vector{Float64}}(:E => zeros(5), :Tg => zeros(5))
            )
            # a ModelDiff carries no model equations of its own
            @test_throws ArgumentError IF.EBMProblem{Float64}(
                IF.ModelDiff{MIZModel,ClassicModel}(), ST, Forcing(0.0), miz_parameters()
            )
        end # @testset "invalid input"

        @testset "show" begin
            prob = EBMProblem(WIModel(); st=(; nx=10, nt=50, dur=1), spectrum=TEST_SPECTRUM)
            @test repr(prob) == "EBMProblem{Float64}(WIModel, SpaceTime{sin}(10, 50, 1), Forcing(0.0) (constant forcing))"
            pretty = sprint(show, MIME"text/plain"(), prob)
            @test occursin("model:      WIModel", pretty)
            @test occursin("parameters: 34 entries", pretty)
            @test occursin("spectrum:   Spectrum(10 components)", pretty)
            @test occursin(
                "spectrum:   nothing",
                sprint(show, MIME"text/plain"(), EBMProblem(ClassicModel(); st=(; nx=5)))
            )
        end # @testset "show"
    end # @testset "EBMProblem"

    @testset "get_diffop" begin
        st = SpaceTime{sin}(60, 100, 1)
        diffop = IF.get_diffop(st)
        @test diffop isa AbstractMatrix{Float64}
        @test size(diffop) == (60, 60)
        # ∂ₓ[(1-x²)∂ₓT] annihilates a constant field
        @test maximum(abs, diffop * ones(60)) < 1e-10
        # and reproduces 2 - 6x² for T = x², except in the polar gridbox where the
        # ghost point closes the no-flux boundary instead
        @test (diffop * (st.x .^ 2))[1:(end - 1)] ≈ (2 .- 6 .* st.x .^ 2)[1:(end - 1)] atol = 1e-8
        # the operator is cached per grid
        @test IF.get_diffop(st) === diffop

        sti = SpaceTime{identity}(40, 100, 1)
        diffopi = IF.get_diffop(sti)
        @test size(diffopi) == (40, 40)
        @test maximum(abs, diffopi * ones(40)) < 1e-10
        @test IF.get_diffop(sti) === diffopi
    end # @testset "get_diffop"

    @testset "hemispheric_mean" begin
        # x is the sine of latitude, so a uniform field has hemispheric mean equal to itself
        x = SpaceTime{sin}(180, 10, 1).x
        @test hemispheric_mean(ones(180), x) ≈ 1.0 rtol = 0.01
        @test hemispheric_mean(fill(-3.5, 180), x) ≈ -3.5 rtol = 0.01
        @test hemispheric_mean(zeros(180), x) == 0.0
        # the mean is linear in the field
        a = @. 7.5 + 20(1 - 2x^2)
        b = @. 3x
        @test hemispheric_mean(a .+ b, x) ≈ hemispheric_mean(a, x) + hemispheric_mean(b, x)
        @test hemispheric_mean(2a, x) ≈ 2hemispheric_mean(a, x)
        # a colder pole lowers the mean
        @test hemispheric_mean(a, x) < hemispheric_mean(fill(a[1], 180), x)
    end # @testset "hemispheric_mean"

    @testset "ice_area" begin
        x = SpaceTime{sin}(180, 10, 1).x
        # full cover integrates to the area of a hemisphere in these units, π/2
        @test ice_area(ones(180), x) ≈ pi / 2 rtol = 0.01
        @test ice_area(zeros(180), x) == 0.0
        # only the polar half is covered
        polar = Float64.(x .> sin(pi / 4))
        @test 0 < ice_area(polar, x) < ice_area(ones(180), x)
        @test ice_area(0.5ones(180), x) ≈ 0.5ice_area(ones(180), x)
        # the documented example
        xdoc = sin.(range(0, pi/2, 181))[1:end-1]
        @test ice_area(@.(2xdoc - xdoc^2), xdoc) ≈ 1.4075808945373096
    end # @testset "ice_area"

    @testset "Solutions storage" begin
        st = SpaceTime{sin}(20, 100, 5)
        par = classic_parameters()
        init = classic_initconds(st)
        vars, sols, annusol = IF.create_storages(
            ClassicModel(), Set((:E, :T)), st, Forcing(0.0), par, init;
            solver=GhostLayerSolver(), lastonly=true
        ) # IF.create_storages
        @test sols isa Solutions{ClassicModel,sin,false}
        @test vars.E == init.E
        @test vars.E !== init.E # create_storages works on a copy of the initial conditions
        @test keys(sols.raw) == Set((:E, :T))
        @test sols.solver isa GhostLayerSolver
        @test sols.lastonly
        @test sols.parameters === par
        @test sols.initconds === init
        # only the final year is stored at full resolution
        @test length(sols.ts) == st.nt
        @test sols.ts[1] ≈ 4 + st.dt / 2
        @test sols.ts[end] ≈ 5 - st.dt / 2
        @test keys(sols.annual) == (:winter, :summer, :avg)
        @test length(sols.annual.avg.E) == st.dur
        # `annusol` is a one-year scratch space for the annual means
        @test annusol.lastonly
        @test length(annusol.ts) == st.nt

        _, fullsols, _ = IF.create_storages(
            ClassicModel(), Set((:E,)), st, Forcing(0.0), par, init;
            solver=NonlinearSolver(), lastonly=false
        ) # IF.create_storages
        @test !fullsols.lastonly
        @test length(fullsols.ts) == st.nt * st.dur
        @test fullsols.ts[1] ≈ st.dt / 2

        @testset "savesol!" begin
            # feed a ramp so that every stored value identifies its own timestep
            for tinx in 1:(st.nt * st.dur)
                IF.savesol!(
                    sols, annusol, Collection{Vector{Float64}}(:E => fill(Float64(tinx), st.nx), :T => zeros(st.nx)),
                    tinx
                ) # IF.savesol!
            end # for tinx
            # raw storage holds the last year only
            @test sols.raw.E[1][1] == 4st.nt + 1
            @test sols.raw.E[end][1] == 5st.nt
            # seasonal peaks are picked out by their within-year index
            @test sols.annual.winter.E[1][1] == st.winter.inx
            @test sols.annual.summer.E[2][1] == st.nt + st.summer.inx
            # and the annual average is the mean over the year
            @test sols.annual.avg.E[1][1] ≈ (1 + st.nt) / 2
            @test sols.annual.avg.E[5][1] ≈ (4st.nt + 1 + 5st.nt) / 2
            # the stored vectors are copies, immune to later mutation of `vars`
            @test sols.raw.E[1] !== sols.raw.E[2]
        end # @testset "savesol!"

        @testset "annual_mean" begin
            means = IF.annual_mean(annusol)
            @test means isa Collection{Vector{Float64}}
            @test keys(means) == Set((:E, :T))
            @test means.E[1] ≈ (4st.nt + 1 + 5st.nt) / 2
            @test means.T == zeros(st.nx)
        end # @testset "annual_mean"

        @testset "savesol! stores copies" begin
            # the solution must not alias the working variables, which are mutated in place
            working = Collection{Vector{Float64}}(:E => fill(1.0, st.nx), :T => zeros(st.nx))
            IF.savesol!(sols, annusol, working, st.nt * st.dur)
            working.E[1] = -999.0
            @test sols.raw.E[end][1] == 1.0
            @test annusol.raw.E[end][1] == 1.0
        end # @testset "savesol! stores copies"
    end # @testset "Solutions storage"

end # @testset "Infrastructure"
