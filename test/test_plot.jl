# Tests for `EnergyBalanceModel.Plot`. No Makie backend is loaded, so the figures are built
# but never rendered; that is enough to exercise the whole plotting path.

@testset "Plot" begin

    @testset "Layout" begin
        layout = Layout([:E :T :h])
        @test layout isa Layout{Symbol}
        @test layout.vars == [:E :T :h]
        @test layout.titles == AbstractString["E" "T" "h"] # titles default to the variable names
        @test size(layout) == (1, 3)
        @test axes(layout, 2) == Base.OneTo(3)
        @test eachindex(layout) == Base.OneTo(3)
        @test layout[1, 2] == (var=:T, title="T")
        @test layout[3] == (var=:h, title="h")

        titled = Layout([:E :T], AbstractString["Enthalpy" "Temperature"])
        @test titled[1, 1].title == "Enthalpy"
        @test_throws ArgumentError Layout([:E :T], AbstractString["one" "two" "three"])
        @test_throws ArgumentError Layout{Symbol}(reshape([:E, :T], 2, 1), AbstractString["a" "b"])
    end # @testset "Layout"

    @testset "default_layout" begin
        @test default_layout(ClassicModel()).vars == [:E :T :h]
        @test default_layout(MIZModel()).vars == [:Ew :Ei :E; :Tw :Ti :T; :h :D :phi]
        @test default_layout(WIModel()).vars == default_layout(MIZModel()).vars
        # each call hands out an independent copy that callers may mutate
        first = default_layout(MIZModel())
        first.titles[1] = "mutated"
        @test default_layout(MIZModel()).titles[1] != "mutated"

        # a difference of solutions gets Δ-prefixed titles
        diff = default_layout(IF.ModelDiff{MIZModel,WIModel}())
        @test diff.vars == default_layout(MIZModel()).vars
        @test all(occursin.(raw"\Delta", string.(diff.titles)))
        # differences involving the classic model fall back to its smaller layout
        @test default_layout(IF.ModelDiff{MIZModel,ClassicModel}()).vars == [:E :T :h]
        @test default_layout(IF.ModelDiff{ClassicModel,MIZModel}()).vars == [:E :T :h]
    end # @testset "default_layout"

    @testset "backend" begin
        # no backend package is loaded in the test environment
        @test ismissing(backend())
        @test_throws PL.BackendError backend(:GLMakie)
        @test_throws PL.BackendError backend(:CairoMakie)
        @test_throws PL.BackendError backend(:NotABackend)
        @test occursin(
            "No Makie backend is currently loaded", sprint(showerror, PL.BackendError(:missing))
        ) # @test occursin
        @test occursin("import GLMakie", sprint(showerror, PL.BackendError(:GLMakie)))
    end # @testset "backend"

    @testset "matricify" begin
        # a vector of per-timestep fields becomes a time × space matrix
        @test PL.matricify([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]) == [1.0 2.0; 3.0 4.0; 5.0 6.0]
        @test size(PL.matricify([zeros(7) for _ in 1:3])) == (3, 7)
    end # @testset "matricify"

    @testset "limit_size" begin
        xs = collect(range(0.0, 1.0, 11))
        ts = collect(range(0.0, 1.0, 21))
        # everything is kept when both dimensions are within the limits
        @test PL.limit_size(xs, ts) == (collect(1:11), collect(1:21))
        # otherwise the points are thinned uniformly, keeping the end points
        xinx, tinx = PL.limit_size(xs, ts, 3, 4)
        @test length(xinx) == 3
        @test length(tinx) == 4
        @test (first(xinx), last(xinx)) == (1, 11)
        @test (first(tinx), last(tinx)) == (1, 21)
        # ranges select a subset
        @test PL.limit_size(xs, ts, 1000, 1000, (0.2, 0.4), (0.0, 0.5))[1] == [3, 4, 5]

        @test_throws ArgumentError PL.limit_size(xs, ts, 1000, 1000, (2.0, 3.0), (0.0, 1.0))
        @test_throws ArgumentError PL.limit_size(xs, ts, 1000, 1000, (0.0, 1.0), (7.0, 8.0))
        @test_throws ArgumentError PL.limit_size(xs, ts, 1000, 1000, (0.4, 0.2), (0.0, 1.0))
        @test_throws ArgumentError PL.limit_size(xs, ts, 1, 1000)
        @test_throws ArgumentError PL.limit_size(xs, ts, 1000, 0)
        # a range holding a single point cannot be contoured
        @test_logs (:warn, r"Only one space point") PL.limit_size(
            xs, ts, 1000, 1000, (0.0, 0.05), (0.0, 1.0)
        ) # @test_logs
    end # @testset "limit_size"

    @testset "get_levels" begin
        data = [1.0 -2.0; 3.0 NaN]
        # ordinary fields get automatic levels
        @test PL.get_levels(Val(false), Val(false), data)[1] == 21
        # floe size gets its own quadratic level set
        dlevels = PL.get_levels(Val(true), Val(false), data)[1]
        @test issorted(dlevels)
        @test first(dlevels) == 0.0
        @test last(dlevels) == 450
        @test PL.get_levels(Val(true), Val(false), data)[4] == collect(0:100:400)
        # differences are symmetric about zero, and ignore NaNs
        difflevels = PL.get_levels(Val(false), Val(true), data)[1]
        @test difflevels[1] ≈ -3.0
        @test difflevels[end] ≈ 3.0
        @test length(difflevels) == 41
        dlevelsdiff = PL.get_levels(Val(true), Val(true), data)[1]
        @test issorted(dlevelsdiff)
        @test dlevelsdiff[1] == -450
        @test dlevelsdiff[end] == 450
    end # @testset "get_levels"

    @testset "plot_raw" begin
        sols = testsol(:miz_gl)
        @test plot_raw(sols) isa Mk.Figure
        @test plot_raw(testsol(:classic_gl)) isa Mk.Figure
        # a difference of solutions switches to the diverging colour scheme
        @test plot_raw(soldiff(testsol(:miz_gl), testsol(:classic_gl))) isa Mk.Figure
        # plots can be placed into an existing figure
        fig = Mk.Figure()
        @test plot_raw(sols, fig[1, 1]) isa Mk.GridPosition
        # and the layout, the ranges and the resolution are all configurable
        @test plot_raw(sols; layout=Layout([:T :h])) isa Mk.Figure
        @test plot_raw(sols; xsizelim=5, tsizelim=5) isa Mk.Figure
        @test plot_raw(sols; xrange=(0.5, 1.0), trange=extrema(sols.ts)) isa Mk.Figure
        @test plot_raw(sols; inspect=true) isa Mk.Figure
        # ranges outside the stored solution are rejected
        @test_throws ArgumentError plot_raw(sols; trange=(100.0, 200.0))
        @test_throws ArgumentError plot_raw(sols; xrange=(2.0, 3.0))
    end # @testset "plot_raw"

    @testset "plot_avg" begin
        sols = testsol(:miz_gl)
        @test plot_avg(sols) isa Mk.Figure
        @test plot_avg(testsol(:classic_gl)) isa Mk.Figure
        @test plot_avg(soldiff(testsol(:miz_gl), testsol(:classic_gl))) isa Mk.Figure
        fig = Mk.Figure()
        @test plot_avg(sols, fig[1, 1]) isa Mk.GridPosition
        @test plot_avg(sols; layout=Layout([:T :phi])) isa Mk.Figure
        @test plot_avg(sols; trange=(1, sols.spacetime.dur)) isa Mk.Figure
        @test_throws ArgumentError plot_avg(sols; trange=(50, 60))
    end # @testset "plot_avg"

    @testset "plot_seasonal" begin
        sols = testsol(:miz_var)
        @test sols isa Solutions{MIZModel,sin,true}
        @test plot_seasonal(sols) isa Mk.Figure
        fig = Mk.Figure()
        @test plot_seasonal(sols, fig[1, 1]) isa Mk.GridPosition
        @test plot_seasonal(sols; title="Custom", xlabel="x", ylabel="y") isa Mk.Figure
        @test plot_seasonal(sols; inspect=true) isa Mk.Figure
        # the axes are user-definable
        @test plot_seasonal(
            testsol(:classic_var);
            yfunc=(s, season, year) -> hemispheric_mean(getproperty(s.annual, season).E[year], s.spacetime.x),
            xfunc=(s, year) -> Float64(year)
        ) isa Mk.Figure
        # a seasonal plot needs a forcing that varies, otherwise there is nothing to trace
        @test_throws MethodError plot_seasonal(testsol(:miz_gl))
    end # @testset "plot_seasonal"

end # @testset "Plot"
