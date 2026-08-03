# Shared helpers and fixtures for the test suite.
#
# This file is included at the top level of `Main` before any `@testset`, so everything
# defined here is visible to every test file.

using EnergyBalanceModel, Test

import LinearAlgebra as LA
import Makie as Mk

# submodule aliases for testing unexported internals
const EBM = EnergyBalanceModel
const IF = EBM.Infrastructure
const UT = EBM.Utilities
const MZ = EBM.MIZEBM
const CL = EBM.ClassicEBM
const WI = EBM.WIMEBM
const PL = EBM.Plot

# Capture everything written to `stdout` by `f` and return it as a `String`. `IOBuffer`s
# cannot be passed to `redirect_stdout`, hence the temporary file.
function capture_stdout(f::Function)::String
    return mktemp() do _, io
        redirect_stdout(io) do
            f()
        end # do
        flush(io)
        seekstart(io)
        return read(io, String)
    end # do
end # function capture_stdout

# The shared space-time grid for all integration tests. It is much coarser than a
# production run (`SpaceTime{sin}(180, 2000, 50)`) but keeps the seasonal ice cycle and the
# model differences resolved, while running in a fraction of a second per model-year.
const ST = SpaceTime{sin}(45, 200, 3)

# uniform 17°C, ice free
const T17 = fill(17.0, ST.nx)

classic_parameters()::Collection{Float64} = default_parameters(ClassicModel())
miz_parameters()::Collection{Float64} = default_parameters(MIZModel())
wi_parameters()::Collection{Float64} = default_parameters(WIModel())

classic_initconds(st::SpaceTime=ST, par::Collection=classic_parameters())::Collection{Vector{Float64}} =
    Collection{Vector{Float64}}(:E => par.cw * fill(17.0, st.nx), :Tg => fill(17.0, st.nx))

function miz_initconds(st::SpaceTime=ST, par::Collection=miz_parameters())::Collection{Vector{Float64}}
    T = fill(17.0, st.nx)
    return Collection{Vector{Float64}}(
        :Ei => zeros(st.nx),
        :Ew => par.cw * T,
        :h => zeros(st.nx),
        :D => zeros(st.nx),
        :Tg => T
    ) # Collection
end # function miz_initconds

# A reduced incident wave spectrum. `initialise(::WIModel, ...)` tabulates the ice
# wavenumber over 10⁵ thicknesses for every frequency, so the cost of a `WIModel` run is
# dominated by the number of spectral components; 10 keeps the tests quick while exercising
# exactly the same code paths as the 30-component default `bretschneider(3.0, 9.5)`.
const TEST_SPECTRUM = bretschneider(3.0, 9.5, collect(range(2pi/23.8, 2pi/2.5; length=10)))

# A forcing that warms and then cools within the simulated period, used for the tests that
# need `Forcing{true}` solutions: base 0 W m⁻² for 1 y, warming to 1 W m⁻² over 1 y, holding
# for 1 y, then cooling to -1 W m⁻² over 2 y. `domain == (0, 1, 2, 3, 5)`.
const VARFORCING = Forcing(0.0, 1.0, -1.0, (1, 1), (1.0, -1.0))

# Integrations are shared between test files: each one is run at most once, on first use.
const SOLUTIONS = Dict{Symbol,Any}()

"""
    testsol(key::Symbol) -> Solutions

Return the shared solution registered under `key`, integrating it on first use.

| key | model | solver | forcing | notes |
|:----|:------|:-------|:--------|:------|
| `:classic_gl`   | `ClassicModel` | `GhostLayerSolver` | constant | |
| `:classic_nl`   | `ClassicModel` | `NonlinearSolver`  | constant | |
| `:classic_full` | `ClassicModel` | `GhostLayerSolver` | constant | `lastonly=false` |
| `:miz_gl`       | `MIZModel`     | `GhostLayerSolver` | constant | |
| `:miz_nl`       | `MIZModel`     | `NonlinearSolver`  | constant | |
| `:wi_gl`        | `WIModel`      | `GhostLayerSolver` | constant | `TEST_SPECTRUM` |
| `:miz_var`      | `MIZModel`     | `GhostLayerSolver` | `VARFORCING` | 5 y, coarse grid |
| `:classic_var`  | `ClassicModel` | `GhostLayerSolver` | `VARFORCING` | 5 y, coarse grid |
"""
function testsol(key::Symbol)
    return get!(SOLUTIONS, key) do
        if key === :classic_gl
            integrate(
                ClassicModel(), ST, Forcing(0.0), classic_parameters(), classic_initconds();
                solver=GhostLayerSolver(), updatefreq=Inf
            )
        elseif key === :classic_nl
            integrate(
                ClassicModel(), ST, Forcing(0.0), classic_parameters(), classic_initconds();
                solver=NonlinearSolver(), updatefreq=Inf
            )
        elseif key === :classic_full
            integrate(
                ClassicModel(), ST, Forcing(0.0), classic_parameters(), classic_initconds();
                solver=GhostLayerSolver(), updatefreq=Inf, lastonly=false
            )
        elseif key === :miz_gl
            integrate(
                MIZModel(), ST, Forcing(0.0), miz_parameters(), miz_initconds();
                solver=GhostLayerSolver(), updatefreq=Inf
            )
        elseif key === :miz_nl
            integrate(
                MIZModel(), ST, Forcing(0.0), miz_parameters(), miz_initconds();
                solver=NonlinearSolver(), updatefreq=Inf
            )
        elseif key === :wi_gl
            integrate(
                WIModel(), ST, Forcing(0.0), wi_parameters(), miz_initconds();
                solver=GhostLayerSolver(), updatefreq=Inf, spectrum=TEST_SPECTRUM
            )
        elseif key === :miz_var
            st = SpaceTime{sin}(20, 100, 5)
            integrate(
                MIZModel(), st, VARFORCING, miz_parameters(), miz_initconds(st);
                solver=GhostLayerSolver(), updatefreq=Inf
            )
        elseif key === :classic_var
            st = SpaceTime{sin}(20, 100, 5)
            integrate(
                ClassicModel(), st, VARFORCING, classic_parameters(), classic_initconds(st);
                solver=GhostLayerSolver(), updatefreq=Inf
            )
        else # unknown key
            throw(ArgumentError("Unknown test solution $key."))
        end # if ===, elseif*7, else
    end # do
end # function testsol

# hemispheric mean of the annual average of `var` in the last simulated year
lastyear_hemi_mean(sols::Solutions, var::Symbol)::Float64 =
    hemispheric_mean(getproperty(sols.annual.avg, var)[sols.spacetime.dur], sols.spacetime.x)
