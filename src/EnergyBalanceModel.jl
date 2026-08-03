"""
    EnergyBalanceModel

A comprehensive package for solving a classic Energy Balance Model (EBM) (Wagner and
Eisenman, 2015) and an extended EBM with the inclusion of a Marginal Ice Zone (MIZ). Other
utilities for data handling and visualization are also provided.

To get started, choose a model (`ClassicModel`, `MIZModel`, or `WIModel`) and bundle it with
its inputs into an `EBMProblem`. Every input is optional and defaults to a sensible value,
so `EBMProblem(model)` already gives a fully-specified, runnable problem; override the
space-time domain, forcing, parameters, initial conditions, or wave spectrum through the
keyword arguments when needed. Then call `solve` to run it. See the documentation of
`EBMProblem`, `solve`, `SpaceTime`, `Forcing`, `default_parameters`, and `Collection` for
details.

The following example runs the wave-influenced EBM (`WIModel`) with all defaults: 50 years
on a 180-point latitudinal grid equally spaced in latitude, 2000 timesteps per year, a
constant forcing of 0.0, initial conditions of uniform temperature 17°C with no ice, and a
default Bretschneider incident wave spectrum. The results are saved and plotted.

```julia-repl
julia> using EnergyBalanceModel

julia> problem = EBMProblem(WIModel())
EBMProblem{Float64} with:
  model:      WIModel
  spacetime:  SpaceTime{sin}(180, 2000, 50)
  forcing:    Forcing(0.0) (constant forcing)
  parameters: 32 entries
  initconds:  5 variables: Set([:Tg, :Ei, :D, :h, :Ew])
  spectrum:   Spectrum(30 components)

julia> sols = solve(problem)
Wavenumber cached in 7.29 s

Integrating WIModel
 100000/100000 [━━━━━━━━━━━━━━━━━━━━━━━━━━━]  100%
 0:39/-0:00 2575.71/sec                     Done ✓
 t = 50.0
Solutions{WIModel, sin, false} with:
  12 solution variables: Set([:Ti, :n, :D, :h, :lambda, :phi, :Ew, :E, :Tw, :T, :Ei, :Ewave])
  on 180 latitudinal gridboxes: [0.00436331, 0.0130896 … 2, 0.999914, 0.99999]
  and 2000 timesteps: 49.00025:0.0005:49.99975
  with forcing Forcing(0.0) (constant forcing)

julia> save(sols, "./example.jld2");

julia> import GLMakie; plot_raw(sols)
```

You can also run a standard example by calling `run_example()`, which uses the `MIZModel` by
default. To run the wave-influenced or classic EBM instead, call `run_example(WIModel())` or
`run_example(ClassicModel())`.

See the documentation for `save`, `load`, `plot_raw`, `plot_avg`, and `plot_seasonal` for
details on data handling and visualisation.
"""
module EnergyBalanceModel

export ClassicModel, MIZModel, WIModel
export Collection, EBMProblem, Forcing, Solutions, SpaceTime
export Spectrum, bretschneider, monochromatic
export default_parameters, integrate, solve
export hemispheric_mean, ice_area
export Layout, backend, default_layout, plot_avg, plot_raw, plot_seasonal
export run_example

include("utilities.jl")
include("infrastructure.jl")
include("mizebm.jl")
include("wimebm.jl")
include("classicebm.jl")
include("plot.jl")

using .ClassicEBM, .Infrastructure, .MIZEBM, .Plot, .Utilities, .WIMEBM

"""
    run_example(model<:AbstractModel=MIZModel(); plotbackend::Symbol=:GLMakie) -> Solutions{M,sin,false}

Run a standard example simulation for the specified `model` (either an instance of
`MIZModel`, `WIModel`, or `ClassicModel`). The results of the last year (year 50) are
plotted using the specified Makie backend `plotbackend`. The backend package must be loaded
beforehand (e.g., `import GLMakie`).

The model is run on a 180-point latitudinal grid equally spaced in latitude, with 2000
timesteps per year for 50 years and a constant forcing of 0.0. The initial conditions are
uniform temperature of 17°C with no ice.

# Examples
```julia-repl
julia> using EnergyBalanceModel, GLMakie

julia> run_example()
Integrating
 100000/100000 [━━━━━━━━━━━━━━━━━━━━━━━━━━━]  100%
 0:15/-0:00 6456.44/sec                     Done ✓
 t = 50.0
Solutions{MIZModel, sin, false} with:
  10 solution variables: Set([:T, :Ei, :Ti, :D, :n, :h, :phi, :Ew, :E, :Tw])
  on 180 latitudinal gridboxes: [0.00436331, 0.0130896 … 2, 0.999914, 0.99999]
  and 2000 timesteps: 49.00025:0.0005:49.99975
  with forcing Forcing{false}(0.0) (constant forcing)

julia> run_example(ClassicModel())
Integrating
 100000/100000 [━━━━━━━━━━━━━━━━━━━━━━━━━━━]  100%
 0:18/-0:00 5702.05/sec                     Done ✓
 t = 50.0
Solutions{ClassicModel, sin, false} with:
  3 solution variables: Set([:T, :h, :E])
  on 180 latitudinal gridboxes: [0.00436331, 0.0130896 … 2, 0.999914, 0.99999]
  and 2000 timesteps: 49.00025:0.0005:49.99975
  with forcing Forcing{false}(0.0) (constant forcing)
```
"""
function run_example(model::M=MIZModel())::Solutions{M,sin,false} where M<:AbstractModel
    problem = EBMProblem(model)
    sols = solve(problem)
    try # plot results
        fig = plot_raw(sols)
        display(fig)
    catch err
        if err isa Plot.BackendError
            msgbuffer = IOBuffer()
            showerror(msgbuffer, err)
            write(msgbuffer, "Skipping plotting.")
            @warn String(take!(msgbuffer))
        else # other error
            rethrow(err)
        end # if isa; else
    end # try; catch
    return sols
end # function run_example

end # module EnergyBalanceModel
