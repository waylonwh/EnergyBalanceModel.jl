# Tests for `EnergyBalanceModel.MIZEBM`, the EBM extended with a marginal ice zone.

@testset "MIZEBM" begin

    st = SpaceTime{sin}(30, 200, 2)
    par = miz_parameters()

    @testset "solar" begin
        x = st.x
        @test MZ.solar(x, 0.0, :ice, par) ≈ @.(par.ai * (par.S0 - par.S1 * x - par.S2 * x^2))
        @test MZ.solar(x, 0.5, :ice, par) ≈ @.(par.ai * (par.S0 + par.S1 * x - par.S2 * x^2))
        @test MZ.solar(x, 0.25, :water, par) ≈ @.((par.a0 - par.a2 * x^2) * (par.S0 - par.S2 * x^2))
        # the Symbol and Val interfaces agree
        @test MZ.solar(x, 0.3, :water, par) == MZ.solar(x, 0.3, Val(:water), par)
        @test MZ.solar(x, 0.3, :ice, par) == MZ.solar(x, 0.3, Val(:ice), par)
        # ice reflects more than open water, and the insolation cycle is annual
        @test all(MZ.solar(x, 0.4, :ice, par) .< MZ.solar(x, 0.4, :water, par))
        @test MZ.solar(x, 0.4, :ice, par) ≈ MZ.solar(x, 1.4, :ice, par)
        # the pole gets its most sunlight in the middle of the year and its least at the ends
        polar = [only(MZ.solar(x[end:end], t, :water, par)) for t in st.t]
        @test abs(st.t[argmax(polar)] - 0.5) <= st.dt
        @test argmin(polar) == 1
        # the equator is far less seasonal than the pole
        equatorial = [only(MZ.solar(x[1:1], t, :water, par)) for t in st.t]
        @test maximum(equatorial) - minimum(equatorial) < maximum(polar) - minimum(polar)
    end # @testset "solar"

    @testset "weighted_avg" begin
        @test MZ.weighted_avg([1.0, 2.0], [3.0, 4.0], [0.0, 1.0]) == [3.0, 2.0]
        @test MZ.weighted_avg([1.0], [3.0], [0.5]) == [2.0]
        @test MZ.weighted_avg(fill(5.0, 3), fill(5.0, 3), [0.0, 0.5, 1.0]) == fill(5.0, 3)
    end # @testset "weighted_avg"

    @testset "temperatures" begin
        @test MZ.water_temp([9.8, -9.8, 0.0], par) == [1.0, -1.0, 0.0]
        @test MZ.water_temp(zeros(3), par) == fill(par.Tm, 3)
        # the ice surface can cool below freezing but never rise above it
        @test MZ.ice_temp([-3.0, 5.0, 0.0], par) == [-3.0, 0.0, 0.0]
    end # @testset "temperatures"

    @testset "wlat" begin
        @test MZ.wlat([par.Tm], par) == [0.0]
        @test MZ.wlat([1.0], par) ≈ [par.m1]
        @test MZ.wlat([2.0], par) ≈ [par.m1 * 2^par.m2]
        # melting is faster in warmer water
        @test only(MZ.wlat([2.0], par)) > only(MZ.wlat([1.0], par))
        # the empirical law is undefined below the freezing point
        @test_throws DomainError MZ.wlat([-1.0], par)
    end # @testset "wlat"

    @testset "concentration" begin
        # φ = -Eᵢ / (Lf h), capped at complete cover
        phi = MZ.concentration([-4.75, -9.5, -19.0, 0.0], [1.0, 1.0, 1.0, 0.0], par)
        @test phi == [0.5, 1.0, 1.0, 0.0]
        @test all(p -> 0 <= p <= 1, phi)
        # gridboxes without ice have zero concentration whatever the enthalpy
        @test MZ.concentration([-5.0, 0.0], zeros(2), par) == zeros(2)
        @test MZ.concentration(zeros(3), ones(3), par) == zeros(3)
        # thicker ice at the same enthalpy means less area covered
        @test only(MZ.concentration([-4.75], [2.0], par)) < only(MZ.concentration([-4.75], [1.0], par))
    end # @testset "concentration"

    @testset "num" begin
        D = [10.0, 10.0, 0.0, 5.0]
        phi = [0.5, 1.0, 0.5, 0.5]
        n = MZ.num(D, phi, par)
        @test n ≈ [0.5, 1.0, 0.0, 0.5] ./ (par.alpha .* [10.0, 10.0, 1.0, 5.0] .^ 2)
        @test n[3] == 0.0 # no floes where there is no floe size
        # smaller floes at the same concentration means more of them
        @test n[4] > n[1]
        @test all(>=(0), n)
    end # @testset "num"

    @testset "area_lead" begin
        D = [10.0, 10.0]
        phi = [0.5, 0.99]
        n = MZ.num(D, phi, par)
        lead = MZ.area_lead(D, phi, n, par)
        @test lead[1] ≈ par.alpha * n[1] * ((D[1] + 2par.rl)^2 - D[1]^2)
        # the lead region can never exceed the open water fraction
        @test lead[2] == 1 - phi[2]
        @test all(lead .<= 1 .- phi)
        @test MZ.area_lead(zeros(2), zeros(2), zeros(2), par) == zeros(2)
    end # @testset "area_lead"

    @testset "fluxes" begin
        @test MZ.bot_flux([1.0, 0.0, 2.0], par) ≈
            par.rhow * par.cp * par.ch * par.u0 .* [1.0, 0.0, 2.0]
        @test MZ.bot_flux([par.Tm], par) == [0.0]

        h = [1.0, 0.0, 1.0]
        D = [10.0, 10.0, 0.0]
        Tw = fill(1.0, 3)
        phi = [1.0, 1.0, 1.0]
        Flat = MZ.lat_flux(h, D, Tw, phi, par)
        @test Flat[1] ≈ phi[1] * h[1] * par.Lf * par.m1 * pi / (par.alpha * D[1])
        @test Flat[2] == 0.0 # no ice thickness, no lateral melt
        @test Flat[3] == 0.0 # no floes, no lateral melt
        @test all(>=(0), Flat)

        # the ghost layer flux relaxes the surface towards Tg, the nonlinear one diffuses
        Tbar = fill(-1.0, st.nx)
        Tg = fill(-2.0, st.nx)
        Fgl = MZ.vert_flux(GhostLayerSolver(), 0.3, :ice, Tbar, 0.0, st, par, Tg)
        @test Fgl ≈ MZ.solar(st.x, 0.3, :ice, par) .- (par.A .+ par.B .* (Tbar .- par.Tm)) .+
            par.cg / par.tau .* (Tg .- Tbar) .+ par.Fb
        Fnl = MZ.vert_flux(NonlinearSolver(), 0.3, :ice, Tbar, 0.0, st, par, Tg)
        @test Fnl ≈ MZ.solar(st.x, 0.3, :ice, par) .- (par.A .+ par.B .* (Tbar .- par.Tm)) .+
            par.D * IF.get_diffop(st) * Tbar .+ par.Fb
        # a uniform surface temperature has no diffusive flux convergence
        @test Fnl ≈ MZ.solar(st.x, 0.3, :ice, par) .- (par.A .+ par.B .* (Tbar .- par.Tm)) .+ par.Fb atol = 1e-8
        # the climate forcing enters additively
        @test MZ.vert_flux(NonlinearSolver(), 0.3, :ice, Tbar, 5.0, st, par, Tg) ≈ Fnl .+ 5.0
    end # @testset "fluxes"

    @testset "redistributeE" begin
        rEi = [-5.0, 3.0, -2.0, 0.0]
        rEw = [4.0, -6.0, 1.0, 0.0]
        Ei, Ew, psiEi, psiEw = MZ.redistributeE(rEi, rEw)
        # ice enthalpy stays non-positive, water enthalpy non-negative
        @test all(<=(0), Ei)
        @test all(>=(0), Ew)
        # and the total enthalpy is conserved
        @test Ei .+ Ew == rEi .+ rEw
        @test Ei == [-5.0, -6.0, -2.0, 0.0]
        @test Ew == [4.0, 3.0, 1.0, 0.0]
        # the transfers are the positive part of Eᵢ and the negative part of E_w
        @test psiEi == [0.0, 3.0, 0.0, 0.0]
        @test psiEw == [0.0, -6.0, 0.0, 0.0]
        @test all(>=(0), psiEi)
        @test all(<=(0), psiEw)
    end # @testset "redistributeE"

    @testset "split_psiEw" begin
        Ql, Qp = MZ.split_psiEw([-2.0, -4.0], [0.5, 1.0], [0.25, 0.5])
        # freezing in the lead region is proportional to its share of the open water
        @test Ql ≈ [-1.0, 0.0]
        @test Qp ≈ [-1.0, -4.0]
        # with complete cover everything is left to pancake growth
        @test Ql .+ Qp == [-2.0, -4.0]
    end # @testset "split_psiEw"

    @testset "dphip and average" begin
        @test MZ.dphip([1.0], par) ≈ [-1 / (par.Lf * par.hmin)]
        @test MZ.dphip([0.0], par) == [0.0]
        # new pancakes of size `fn` are averaged into the existing floes
        @test MZ.average([2.0, 5.0], 1.0, [1.0, 0.0], [1.0, 0.0]) == [1.5, 0.0]
        @test MZ.average([2.0], 1.0, [3.0], [1.0]) ≈ [(3 * 2 + 1 * 1) / 4]
        @test MZ.average([2.0], 1.0, [1.0], [0.0]) == [2.0] # nothing new to average in
    end # @testset "dphip and average"

    @testset "tendencies" begin
        @test MZ.Ei_t([0.5], [2.0], [1.0], [3.0]) == [0.5 * 2 + 1 + 0.5 * 3]
        @test MZ.Ew_t([0.5], [2.0], [1.0], [3.0]) == [0.5 * 2 - 1 - 0.5 * 3]
        # what leaves the water through the lateral and basal fluxes enters the ice
        @test only(MZ.Ei_t([0.5], [0.0], [1.0], [3.0])) == -only(MZ.Ew_t([0.5], [0.0], [1.0], [3.0]))
        @test MZ.h_t([1.0], [1.0], par) ≈ [-2 / par.Lf]
        @test MZ.h_t([0.0], [0.0], par) == [0.0]
        @test MZ.forward_euler([1.0, 2.0], [2.0, -4.0], 0.5) == [2.0, 0.0]

        # D_t = lateral melt (negative) + growth into the leads + welding
        Dt = MZ.D_t(
            [1.0, 1.0], [10.0, 10.0], [-1.0, -1.0], [1.0, 1.0], [0.5, 0.5], [0.0, 0.0], par;
            breakup=BitArray([false, true])
        ) # MZ.D_t
        melt = -pi / 2par.alpha * par.m1
        weld = par.kappa * par.alpha / 4 * 0.5 * 10.0^3
        @test Dt ≈ [melt + weld, melt] # welding is suppressed where the floes are breaking
        # welding also stops once the ice surface reaches the melting point
        @test MZ.D_t(
            [1.0], [10.0], [par.Tm], [1.0], [0.5], [0.0], par; breakup=BitArray([false])
        ) ≈ [melt]
        # freezing in the leads (Ql < 0) grows the floes
        @test only(MZ.D_t(
            [1.0], [10.0], [-1.0], [1.0], [0.5], [-1.0], par; breakup=BitArray([false])
        )) > melt + weld
    end # @testset "tendencies"

    @testset "T0eq and solveT0" begin
        diffop = IF.get_diffop(st)
        h = fill(1.0, st.nx)
        Tw = fill(-1.0, st.nx)
        phi = fill(0.5, st.nx)
        T0 = MZ.solveT0(NonlinearSolver(), zeros(st.nx), st.x, 0.3, h, Tw, phi, diffop, 0.0, par)
        @test length(T0) == st.nx
        @test all(isfinite, T0)
        @test maximum(abs, MZ.T0eq(T0, st.x, 0.3, h, Tw, phi, 0.0, diffop, par)) < 1e-8
        # a NaN guess is replaced by the melting point before the solve
        @test MZ.solveT0(
            NonlinearSolver(), fill(NaN, st.nx), st.x, 0.3, h, Tw, phi, diffop, 0.0, par
        ) ≈ T0
        # ice-free gridboxes (h = 0) fall back to `hmin`, so the solve stays finite
        @test all(isfinite, MZ.solveT0(
            NonlinearSolver(), zeros(st.nx), st.x, 0.3, zeros(st.nx), Tw, phi, diffop, 0.0, par
        ))
        # a stronger forcing warms the ice surface
        T0warm = MZ.solveT0(NonlinearSolver(), zeros(st.nx), st.x, 0.3, h, Tw, phi, diffop, 20.0, par)
        @test all(T0warm .> T0)

        # the ghost layer scheme solves the same balance diagnostically
        T0gl = MZ.solveT0(GhostLayerSolver(), st.x, 0.3, h, fill(-2.0, st.nx), Tw, phi, 0.0, par)
        @test length(T0gl) == st.nx
        @test all(isfinite, T0gl)
    end # @testset "T0eq and solveT0"

    @testset "initialise" begin
        init = miz_initconds(st)
        vars, sols, annusol = IF.initialise(
            MIZModel(), st, Forcing(0.0), miz_parameters(), init;
            solver=GhostLayerSolver(), lastonly=true, spectrum=nothing
        ) # IF.initialise
        @test sols isa Solutions{MIZModel,sin,false}
        @test annusol isa Solutions{MIZModel,sin,false}
        @test keys(sols.raw) == Set((:Ei, :Ew, :D, :h, :E, :Ti, :Tw, :T, :phi, :n))
        # the diagnostics needed by the first `step!` are filled in
        @test keys(vars) ⊇ Set((:phi, :Tw, :T0, :Ti, :T))
        @test vars.phi == zeros(st.nx) # the initial condition carries no ice
        @test vars.Tw ≈ init.Ew ./ miz_parameters().cw
        @test all(isfinite, vars.Ti) # NaNs are removed at initialisation
        @test vars.T ≈ vars.Tw # with no ice the surface is the mixed layer
    end # @testset "initialise"

    @testset "step!" for solver in (GhostLayerSolver(), NonlinearSolver())
        stepst = SpaceTime{sin}(30, 200, 3)
        vars, _, _ = IF.initialise(
            MIZModel(), stepst, Forcing(0.0), par, miz_initconds(stepst);
            solver, lastonly=true, spectrum=nothing
        ) # IF.initialise
        @test IF.step!(
            MIZModel(), stepst.t[1], 0.0, vars, stepst, par; solver, spectrum=nothing
        ) === vars

        # run three years so that a seasonal ice cover develops
        for ti in 2:(stepst.nt * stepst.dur)
            IF.step!(
                MIZModel(), stepst.t[mod1(ti, stepst.nt)], 0.0, vars, stepst, par;
                solver, spectrum=nothing
            ) # IF.step!
        end # for ti

        @test any(>(0), vars.phi) # ice forms
        @test all(isfinite, vars.Ei)
        @test all(isfinite, vars.Ew)
        @test all(isfinite, vars.h)
        @test all(isfinite, vars.D)
        @test all(isfinite, vars.T)
        # the enthalpy split is exhaustive and correctly signed
        @test vars.E ≈ vars.Ei .+ vars.Ew
        @test all(<=(0), vars.Ei)
        @test all(>=(0), vars.Ew)
        # the sea ice state stays inside its physical bounds
        @test all(p -> 0 <= p <= 1, vars.phi)
        @test all(>=(0), vars.h)
        @test all(d -> 0 <= d <= par.Dmax, vars.D)
        # ice thickness, floe size and concentration all vanish together
        @test all(iszero, vars.h[iszero.(vars.Ei)])
        @test all(iszero, vars.D[iszero.(vars.Ei)])
        @test all(iszero, vars.phi[iszero.(vars.Ei)])
        # the ice surface temperature is NaN exactly where there is no ice
        @test all(isnan, vars.Ti[iszero.(vars.Ei)])
        @test all(isfinite, vars.Ti[vars.Ei .< 0])
        @test all(vars.Ti[vars.Ei .< 0] .<= par.Tm)
        # the diagnostics are consistent with the prognostic state
        @test vars.Tw ≈ MZ.water_temp(vars.Ew, par)
        @test vars.phi ≈ MZ.concentration(vars.Ei, vars.h, par)
        @test vars.n ≈ MZ.num(vars.D, vars.phi, par)
        @test vars.T ≈ MZ.weighted_avg(replace(vars.Ti, NaN => 0.0), vars.Tw, vars.phi)
    end # @testset "step!"

    @testset "breakup keyword" begin
        stepst = SpaceTime{sin}(30, 200, 3)
        vars, _, _ = IF.initialise(
            MIZModel(), stepst, Forcing(0.0), par, miz_initconds(stepst);
            solver=GhostLayerSolver(), lastonly=true, spectrum=nothing
        ) # IF.initialise
        for ti in 1:(stepst.nt * stepst.dur)
            IF.step!(
                MIZModel(), stepst.t[mod1(ti, stepst.nt)], 0.0, vars, stepst, par;
                solver=GhostLayerSolver(), spectrum=nothing
            ) # IF.step!
        end # for ti
        # welding is suppressed in the gridboxes flagged as breaking, so the floes stay smaller
        welded = deepcopy(vars)
        broken = deepcopy(vars)
        IF.step!(
            MIZModel(), stepst.t[1], 0.0, welded, stepst, par;
            solver=GhostLayerSolver(), spectrum=nothing
        ) # IF.step!
        IF.step!(
            MIZModel(), stepst.t[1], 0.0, broken, stepst, par;
            solver=GhostLayerSolver(), spectrum=nothing, breakup=trues(stepst.nx)
        ) # IF.step!
        @test all(broken.D .<= welded.D)
        @test any(broken.D .< welded.D)
    end # @testset "breakup keyword"

    @testset "unimplemented solver" begin
        # `ActiveSetSolver` has no MIZEBM implementation yet; see test_solvers.jl
        @test_throws ErrorException MZ.step_temperature!(
            ActiveSetSolver(), Collection{Vector{Float64}}()
        ) # @test_throws
    end # @testset "unimplemented solver"

end # @testset "MIZEBM"
