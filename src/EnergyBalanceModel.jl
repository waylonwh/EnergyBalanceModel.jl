"""
    EnergyBalanceModel

A comprehensive package for solving a classic Energy Balance Model (EBM) (Wagner and
Eisenman, 2015) and an extended EBM with the inclusion of a Marginal Ice Zone (MIZ). Other
utilities for data handling and visualization are also provided.

To get started, define a space-time domain, a forcing function, parameters, and initial
conditions. Then call `integrate` to run the model. See documentation of `SpaceTime`,
`Forcing`, `default_parameters`, `Collection`, and `integrate` for details. The following
example runs the EBM with MIZ for 30 years on a 180-point latitudinal grid equally spaced in
sine latitude, with 2000 timesteps per year and a constant forcing of 0.0. Then the results
are saved and plotted.

```julia-repl
julia> using EnergyBalanceModel

julia> st = SpaceTime{sin}(180, 2000, 30)
SpaceTime{sin} with:
  180 latitudinal gridboxes: [0.00436331, 0.0130896, … 762, 0.999914, 0.99999]
  2000 timesteps per year: [0.00025, 0.00075, 0.001 … 99875, 0.99925, 0.99975]
  30 years of simulation: t∈[0,30]
  winter at t=0.26125, summer at t=0.77375

julia> forcing = Forcing(0.0)
Forcing{true}(0.0) is constant:
  F(t)=0.0, t∈[0,∞)

julia> par = default_parameters(MIZ)
Collection{Float64} with 24 entries:
  :Dmax  => 156.0
  :a2    => 0.1
  :alpha => 0.66
  :m1    => 50.4576
  :D     => 0.6
  :S1    => 338.0
  :B     => 2.1
  ⋮      => ⋮

julia> T = fill(17.0, st.nx);

julia> init = Collection{Vec}(
           :Ei => zeros(st.nx),
           :Ew => T .* par.cw,
           :h => zeros(st.nx),
           :D => zeros(st.nx),
           :Tg => T,
       )
Collection{Vector{Float64}} with 5 entries:
  :Tg => [17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0  …  17.0, …
  :Ei => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0…
  :D  => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0…
  :h  => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0…
  :Ew => [166.6, 166.6, 166.6, 166.6, 166.6, 166.6, 166.6, 166.6, 166.6, 166.6 …

julia> sols = integrate(MIZ, st, forcing, par, init)
Integrating
 60000/60000 [━━━━━━━━━━━━━━━━━━━━━━━━━━━━━]  100%
 0:46/-0:00 1295.08/sec                     Done ✓
 t = 30.0
Solutions{EnergyBalanceModel.Infrastructure.MIZModel, sin, true} with:
  10 solution variables: Set([:T, :Ei, :Ti, :D, :n, :h, :phi, :Ew, :E, :Tw])
  on 180 latitudinal gridboxes: [0.00436331, 0.0130896 … 2, 0.999914, 0.99999]
  and 2000 timesteps: 29.00025:0.0005:29.99975
  with forcing Forcing{true}(0.0) (constant forcing)

julia> save(sols, "./miz_sol.jld2")
"./miz_sol.jld2"

julia> import GLMakie; plot_raw(sols)
```

See the documentation for submodules `IO` and `Plot` for details on data handling and
visualization.
"""
module EnergyBalanceModel

export MIZ, Classic
export Vec, Collection, SpaceTime, Forcing, Solutions
export integrate, default_parameters
export annual_mean, hemispheric_mean
export safehouse, house!, retrieve, save, load!
export Layout, backend, plot_raw, plot_avg, plot_seasonal

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

end # module EnergyBalanceModel
