"""
    EnergyBalanceModel

A comprehensive package for solving a classic Energy Balance Model (EBM) (Wagner and
Eisenman, 2015) and an extended EBM with the inclusion of a Marginal Ice Zone (MIZ). Other
utilities for data handling and visualization are also provided.

To get started, define a space-time domain, a forcing function, parameters, and initial
conditions. Then call `integrate` to run the model. See documentation of `SpaceTime`,
`Forcing`, `default_parameters`, `Collection`, and `integrate` for details. The following
example runs the EBM with MIZ for 50 years on a 180-point latitudinal grid equally spaced in
latitude, with 2000 timesteps per year and a constant forcing of 0.0. The initial conditions
are uniform temperature of 17°C with no ice. The results are saved and plotted.

```julia-repl
julia> using EnergyBalanceModel

julia> st = SpaceTime{sin}(180, 2000, 50);

julia> forcing = Forcing(0.0);

julia> par = default_parameters(miz);

julia> T = fill(17.0, st.nx);

julia> init = Collection{Vec}(
           :Ei => zeros(st.nx),
           :Ew => T .* par.cw,
           :h => zeros(st.nx),
           :D => zeros(st.nx),
           :Tg => T,
       );

julia> sols = integrate(miz, st, forcing, par, init)
Integrating
 100000/100000 [━━━━━━━━━━━━━━━━━━━━━━━━━━━]  100%
 1:21/-0:00 1231.41/sec                     Done ✓
 t = 50.0
Solutions{EnergyBalanceModel.Infrastructure.MIZModel, sin, true} with:
  10 solution variables: Set([:T, :Ei, :Ti, :D, :n, :h, :phi, :Ew, :E, :Tw])
  on 180 latitudinal gridboxes: [0.00436331, 0.0130896 … 2, 0.999914, 0.99999]
  and 2000 timesteps: 49.00025:0.0005:49.99975
  with forcing Forcing{true}(0.0) (constant forcing)

julia> save(sols, "./example.jld2");

julia> import GLMakie; plot_raw(sols)
```

You can also run the example above by calling `run_example()`. To run the classic EBM
instead, call `run_example(classic)`.

See the documentation for `save`, `load`, `plot_raw`, `plot_avg`, and `plot_seasonal` for
details on data handling and visualisation.
"""
module EnergyBalanceModel

export miz, classic
export Vec, Collection, SpaceTime, Forcing, Solutions
export integrate, default_parameters
export annual_mean, hemispheric_mean
export safehouse, house!, retrieve, save, load!
export Layout, backend, plot_raw, plot_avg, plot_seasonal
export run_example

include("utilities.jl")
include("infrastructure.jl")
include("mizebm.jl")
include("classicebm.jl")
include("plot.jl")
include("io.jl")

using .Utilities
using .Infrastructure
using .MIZEBM
using .ClassicEBM
using .Plot
using .IO

"""
    run_example(model<:AbstractModel=miz; saveto::String="./example.jld2", plotbackend::Symbol=:GLMakie)

Run a standard example simulation for the specified `model` (either `miz` or `classic`).
The results are saved to the file path given by `saveto`. The results of the last year
(year 50) are plotted using the specified Makie backend `plotbackend`. The backend package
must be loaded beforehand (e.g., `import GLMakie`).

The model is run on a 180-point latitudinal grid equally spaced in sine latitude, with 2000
timesteps per year for 50 years and a constant forcing of 0.0. The initial conditions are
uniform temperature of 17°C with no ice.

# Examples
```julia-repl
julia> using EnergyBalanceModel

julia> import GLMakie; run_example()
Integrating
 100000/100000 [━━━━━━━━━━━━━━━━━━━━━━━━━━━]  100%
 1:15/-0:00 1335.34/sec                     Done ✓
 t = 50.0
Solutions{EnergyBalanceModel.Infrastructure.MIZModel, identity, true} with:
  10 solution variables: Set([:T, :Ei, :Ti, :D, :n, :h, :phi, :Ew, :E, :Tw])
  on 180 latitudinal gridboxes: [0.00277778, 0.0083333 … , 0.991667, 0.997222]
  and 2000 timesteps: 49.00025:0.0005:49.99975
  with forcing Forcing{true}(0.0) (constant forcing)

julia> run_example(classic)
Integrating
 100000/100000 [━━━━━━━━━━━━━━━━━━━━━━━━━━━]  100%
 0:17/-0:00 6008.36/sec                     Done ✓
 t = 50.0
┌ Warning: File ./example.jld2 already exists. Last modified on 7 Nov 2025 at 15:38:39. The EXISTING file has been renamed to ./example_ce5a453c.jld2.
└ @ EnergyBalanceModel.IO EnergyBalanceModel.jl/src/io.jl:48
Solutions{EnergyBalanceModel.Infrastructure.ClassicModel, identity, true} with:
  3 solution variables: Set([:T, :h, :E])
  on 180 latitudinal gridboxes: [0.00277778, 0.0083333 … , 0.991667, 0.997222]
  and 2000 timesteps: 49.00025:0.0005:49.99975
  with forcing Forcing{true}(0.0) (constant forcing)
```
"""
function run_example(
    model::M=miz;
    saveto::String="./example.jld2", plotbackend::Symbol=:GLMakie
) where M<:AbstractModel
    st = SpaceTime(180, 2000, 50)
    forcing = Forcing(0.0)
    par = default_parameters(model)
    T = fill(17.0, st.nx)
    init = Collection{Vec}(:Tg => T)
    if model isa MIZModel
        init.Ei = zeros(st.nx)
        init.Ew = T .* par.cw
        init.h = zeros(st.nx)
        init.D = zeros(st.nx)
    elseif model isa ClassicModel
        init.E = par.cg .* T
    # no else since default_parameters would error earlier
    end # if isa
    sols = integrate(model, st, forcing, par, init)
    save(sols, saveto)
    try # plot results
        plot_raw(sols, plotbackend)
    catch e
        if startswith(e.msg, "Backend not loaded.") # backend error
            @warn string(e.msg, " Skipping plotting of results.")
        else # other error
            rethrow(e)
        end # if startswith; else
    end # try; catch
    return sols
end # function run_example

const VER = 191522

end # module EnergyBalanceModel
