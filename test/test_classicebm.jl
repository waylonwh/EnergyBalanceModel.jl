# Tests for `EnergyBalanceModel.ClassicEBM`, the Wagner & Eisenman (2015) model.

@testset "ClassicEBM" begin

    st = SpaceTime{sin}(30, 200, 2)
    par = classic_parameters()

    @testset "get_statics" begin
        stat = CL.get_statics(st, par)
        @test keys(stat) == (:cg_tau, :dt_tau, :dc, :kappa, :S, :M, :aw, :kLf)
        @test stat.cg_tau ≈ par.cg / par.tau
        @test stat.dt_tau ≈ st.dt / par.tau
        @test stat.dc ≈ stat.dt_tau * stat.cg_tau
        @test stat.M ≈ par.B + par.cg / par.tau
        @test stat.kLf ≈ par.k * par.Lf
        @test stat.aw ≈ @. par.a0 - par.a2 * st.x^2
        @test size(stat.kappa) == (st.nx, st.nx)
        @test stat.kappa ≈ (1 + stat.dt_tau) * LA.I(st.nx) -
            st.dt * par.D * IF.get_diffop(st) / par.cg

        # the insolation matrix carries one extra column wrapping around to the new year
        @test size(stat.S) == (st.nx, st.nt + 1)
        @test stat.S[:, end] == stat.S[:, 1]
        for j in (1, 50, st.nt)
            @test stat.S[:, j] ≈ @. par.S0 - par.S2 * st.x^2 - par.S1 * cos(2pi * st.t[j]) * st.x
        end # for j
        # insolation at the pole peaks in the middle of the year and bottoms out at its ends
        @test abs(st.t[argmax(stat.S[end, 1:st.nt])] - 0.5) <= st.dt
        @test argmin(stat.S[end, 1:st.nt]) == 1

        # results are cached and only recomputed when the grid or the parameters change
        @test CL.get_statics(st, par).kappa === stat.kappa
        @test CL.get_statics(st, par).S === stat.S
        other = classic_parameters()
        other.D = 2par.D
        @test CL.get_statics(st, other).kappa !== stat.kappa
    end # @testset "get_statics"

    @testset "surface_temperature" begin
        # ice-free gridboxes (E >= 0) follow the mixed layer, T = E/cw
        @test CL.surface_temperature([9.8], [-5.0], par) == [1.0]
        @test CL.surface_temperature([0.0], [-5.0], par) == [0.0]
        # ice-covered gridboxes take the surface temperature, but never above freezing
        @test CL.surface_temperature([-1.0], [-2.0], par) == [-2.0]
        @test CL.surface_temperature([-1.0], [3.0], par) == [0.0]
        @test CL.surface_temperature(
            [9.8, -1.0, -1.0], [0.0, -2.0, 3.0], par
        ) == [1.0, -2.0, 0.0]
    end # @testset "surface_temperature"

    @testset "T0eq and solveT0" begin
        diffop = IF.get_diffop(st)
        aS = MZ.solar(st.x, 0.3, :ice, par)
        E = [fill(-5.0, 15); fill(50.0, 15)] # ice at the equator half, open water at the pole half
        h = @. -E / par.Lf * (E < 0)

        # where there is no ice the equation degenerates to T0 = 0
        @test CL.T0eq(fill(2.0, st.nx), E, h, aS, diffop, 0.0, par)[16:end] == fill(2.0, 15)

        T0 = CL.solveT0(zeros(st.nx), E, h, aS, diffop, 0.0, par)
        @test length(T0) == st.nx
        @test maximum(abs, CL.T0eq(T0, E, h, aS, diffop, 0.0, par)) < 1e-8
        @test all(iszero, T0[16:end])
        @test all(isfinite, T0)
        # a stronger forcing warms the ice surface
        T0warm = CL.solveT0(zeros(st.nx), E, h, aS, diffop, 20.0, par)
        @test all(T0warm[1:15] .> T0[1:15])
    end # @testset "T0eq and solveT0"

    @testset "initialise" begin
        init = classic_initconds(st)
        vars, sols, annusol = IF.initialise(
            ClassicModel(), st, Forcing(0.0), par, init;
            solver=GhostLayerSolver(), lastonly=true, spectrum=nothing
        ) # IF.initialise
        @test keys(sols.raw) == Set((:E, :T, :h))
        @test keys(vars) == Set((:E, :Tg))
        @test vars.E == init.E
        @test sols isa Solutions{ClassicModel,sin,false}
        @test annusol isa Solutions{ClassicModel,sin,false}
    end # @testset "initialise"

    @testset "step!" for solver in (GhostLayerSolver(), NonlinearSolver())
        init = classic_initconds(st)
        vars, _, _ = IF.initialise(
            ClassicModel(), st, Forcing(0.0), par, init;
            solver, lastonly=true, spectrum=nothing
        ) # IF.initialise
        stepped = IF.step!(ClassicModel(), st.t[1], 0.0, vars, st, par; solver, spectrum=nothing)
        @test stepped === vars
        @test keys(vars) ⊇ Set((:E, :T, :h))
        @test all(isfinite, vars.E)
        @test all(isfinite, vars.T)
        # the ice thickness is diagnosed from the enthalpy
        @test vars.h == @.(-vars.E / par.Lf * (vars.E < 0))
        @test all(>=(0), vars.h)
        # a warm ice-free start stays ice free after one step
        @test all(iszero, vars.h)

        # ice grows where the enthalpy is negative
        cold = Collection{Vector{Float64}}(:E => fill(-9.5, st.nx), :Tg => fill(-5.0, st.nx))
        coldvars, _, _ = IF.initialise(
            ClassicModel(), st, Forcing(0.0), par, cold; solver, lastonly=true, spectrum=nothing
        ) # IF.initialise
        for ti in 1:st.nt
            IF.step!(ClassicModel(), st.t[ti], 0.0, coldvars, st, par; solver, spectrum=nothing)
        end # for ti
        @test all(isfinite, coldvars.E)
        @test all(isfinite, coldvars.T)
        @test any(>(0), coldvars.h) # a year of seasons does not melt the whole hemisphere
        # `T` is diagnosed from the enthalpy at the start of the step [WE15 Eq. (9)]
        before = copy(coldvars.E)
        IF.step!(ClassicModel(), st.t[1], 0.0, coldvars, st, par; solver, spectrum=nothing)
        @test all(coldvars.T[before .< 0] .<= 0) # ice is never above freezing
        @test coldvars.T[before .>= 0] ≈ before[before .>= 0] ./ par.cw
    end # @testset "step!"

    @testset "unimplemented solver" begin
        # `ActiveSetSolver` has no ClassicEBM implementation yet; see test_solvers.jl
        @test_throws ErrorException CL.step_temperature!(
            ActiveSetSolver(), Collection{Vector{Float64}}()
        ) # @test_throws
    end # @testset "unimplemented solver"

end # @testset "ClassicEBM"
