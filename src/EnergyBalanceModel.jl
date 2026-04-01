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

julia> par = default_parameters(MIZModel());

julia> T = fill(17.0, st.nx);

julia> init = Collection{Vec}(
           :Ei => zeros(st.nx),
           :Ew => T .* par.cw,
           :h => zeros(st.nx),
           :D => zeros(st.nx),
           :Tg => T,
       );

julia> sols = integrate(MIZModel(), st, forcing, par, init)
Integrating
 100000/100000 [━━━━━━━━━━━━━━━━━━━━━━━━━━━]  100%
 1:21/-0:00 1231.41/sec                     Done ✓
 t = 50.0
Solutions{EnergyBalanceModel.Infrastructure.MIZModel, sin, false} with:
  10 solution variables: Set([:T, :Ei, :Ti, :D, :n, :h, :phi, :Ew, :E, :Tw])
  on 180 latitudinal gridboxes: [0.00436331, 0.0130896 … 2, 0.999914, 0.99999]
  and 2000 timesteps: 49.00025:0.0005:49.99975
  with forcing Forcing{false}(0.0) (constant forcing)

julia> save(sols, "./example.jld2");

julia> import GLMakie; plot_raw(sols)
```

You can also run the example above by calling `run_example()`. To run the classic EBM
instead, call `run_example(ClassicModel())`.

See the documentation for `save`, `load`, `plot_raw`, `plot_avg`, and `plot_seasonal` for
details on data handling and visualisation.
"""
module EnergyBalanceModel

export ClassicModel, MIZModel
export Collection, Forcing, Par, Solutions, SpaceTime, Vec
export default_parameters, integrate
export hemispheric_mean, ice_area
export Layout, backend, plot_avg, plot_raw, plot_seasonal
export run_example

include("utilities.jl")
include("infrastructure.jl")
include("mizebm.jl")
include("classicebm.jl")
include("plot.jl")

using .ClassicEBM, .Infrastructure, .MIZEBM, .Plot, .Utilities

"""
    run_example(model<:AbstractModel=MIZModel(); plotbackend::Symbol=:GLMakie) -> Solutions{M,sin,false}

Run a standard example simulation for the specified `model` (either an instance of
`MIZModel` or `ClassicModel`). The results of the last year (year 50) are plotted using the
specified Makie backend `plotbackend`. The backend package must be loaded beforehand
(e.g., `import GLMakie`).

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
function run_example(
    model::M=MIZModel(); plotbackend::Symbol=Plot.find_backend()
)::Solutions{M,sin,false} where M<:AbstractModel
    st = SpaceTime{sin}(180, 2000, 50)
    forcing = Forcing(0.0)
    par = default_parameters(model)
    T = fill(17.0, st.nx)
    init = Collection{Vec}(:Tg => T)
    if model isa MIZModel
        init.Ei = zeros(st.nx)
        init.Ew = par.cw * T
        init.h = zeros(st.nx)
        init.D = zeros(st.nx)
    elseif model isa ClassicModel
        init.E = par.cw * T
    # no else since default_parameters would error earlier
    end # if isa; elseif
    sols = integrate(model, st, forcing, par, init)
    try # plot results
        fig = plot_raw(sols, plotbackend)
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

import PrecompileTools as PT

PT.@setup_workload begin
    import InteractiveUtils as IU
    ms = Tuple(M() for M in IU.subtypes(AbstractModel) if M !== Infrastructure.ModelDiff)
    Fs = (identity, sin)
    fs_args = ((0.0,), (0.0, 1.0, 0.0, (1, 1), (1.0, -1.0)))
    m2s = Dict{AbstractModel,Solutions}()
    redirect_stdout(devnull)
    redirect_stderr(devnull)
    PT.@compile_workload begin
        for m in ms, F in Fs, farg in fs_args
            st = SpaceTime{F}(10, 10, 1)
            forcing = Forcing(farg...)
            par = default_parameters(m)
            T = fill(0.0, st.nx)
            init = Collection{Vec}(:Tg => T)
            if m isa ClassicModel
                init.E = par.cw * T
            else # MIZModel
                init.Ei = zeros(st.nx)
                init.Ew = par.cw * T
                init.h = zeros(st.nx)
                init.D = zeros(st.nx)
            end # if isa; elseif
            sol = integrate(m, st, forcing, par, init)
            F === sin && length(farg) == 1 && (m2s[m] = sol)
        end # for m, F, farg
        m2s[MIZModel()] - m2s[ClassicModel()]
    end # PT.@compile_workload begin
end # PT.@setup_workload begin

end # module EnergyBalanceModel
