# Tests for the solvers that advance the surface temperature.

@testset "Solvers" begin

    @testset "the solver is recorded in the solution" begin
        @test testsol(:classic_gl).solver isa GhostLayerSolver
        @test testsol(:classic_nl).solver isa NonlinearSolver
        @test testsol(:miz_gl).solver isa GhostLayerSolver
        @test testsol(:miz_nl).solver isa NonlinearSolver
        # a difference of two solutions keeps both provenances
        @test soldiff(testsol(:miz_gl), testsol(:miz_nl)).solver isa
            IF.DiffSolver{GhostLayerSolver,NonlinearSolver}
        @test soldiff(testsol(:miz_gl), testsol(:classic_gl)).solver isa GhostLayerSolver
        @test occursin(
            "solved by GhostLayerSolver", sprint(show, MIME"text/plain"(), testsol(:miz_gl))
        ) # @test occursin
    end # @testset "the solver is recorded in the solution"

    @testset "the schemes agree on the climate" begin
        # the ghost layer and the nonlinear scheme discretise the same balance, so their
        # equilibrium climates must be close even though the timestepping differs
        for var in (:T, :E)
            @test abs(
                lastyear_hemi_mean(testsol(:classic_gl), var) -
                lastyear_hemi_mean(testsol(:classic_nl), var)
            ) < 1.0
        end # for var
        @test abs(
            lastyear_hemi_mean(testsol(:miz_gl), :T) - lastyear_hemi_mean(testsol(:miz_nl), :T)
        ) < 1.0
        @test abs(
            lastyear_hemi_mean(testsol(:miz_gl), :E) - lastyear_hemi_mean(testsol(:miz_nl), :E)
        ) < 10.0
    end # @testset "the schemes agree on the climate"

    @testset "ActiveSetSolver is not implemented yet" begin
        # `ActiveSetSolver` is the documented default of `integrate` and `solve`, but no
        # model implements `step_temperature!` for it yet. These tests are expected to fail
        # and will report an error once the scheme lands, as a reminder to un-break them.
        st = SpaceTime{sin}(10, 50, 1)
        @test_broken integrate(
            ClassicModel(), st, Forcing(0.0), classic_parameters(), classic_initconds(st);
            solver=ActiveSetSolver(), updatefreq=Inf
        ) isa Solutions
        @test_broken integrate(
            MIZModel(), st, Forcing(0.0), miz_parameters(), miz_initconds(st);
            solver=ActiveSetSolver(), updatefreq=Inf
        ) isa Solutions
        @test_broken integrate(
            WIModel(), st, Forcing(0.0), wi_parameters(), miz_initconds(st);
            solver=ActiveSetSolver(), updatefreq=Inf, spectrum=TEST_SPECTRUM
        ) isa Solutions
        # ... which is also what the default solver of `solve` currently does
        @test_broken solve(
            EBMProblem(ClassicModel(); st=(; nx=10, nt=50, dur=1)); updatefreq=Inf
        ) isa Solutions
    end # @testset "ActiveSetSolver is not implemented yet"

end # @testset "Solvers"
