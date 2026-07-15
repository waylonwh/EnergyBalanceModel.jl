module Plot

using ..Infrastructure, ..Utilities

import InteractiveUtils as IU, Makie as Mk

export Layout, default_layout, backend
export plot_avg, plot_raw, plot_seasonal

"""
    Layout(vars::Matrix{T}, titles::Matrix{AbstractString})

A layout structure for plotting. When using `plot_raw` or `plot_avg`, type `T` should be
`Symbol`, representing variable names in `Solutions`. The `titles` matrix contains
corresponding titles for each subplot.

# Examples
```julia-repl
julia> Layout([:E :T :h], AbstractString["Enthalpy" "Temperature" "Thinkness"])
Layout{Symbol}([:E :T :h], AbstractString["Enthalpy" "Temperature" "Thinkness"])
```
"""
struct Layout{T}
    vars::Matrix{T}
    titles::Matrix{AbstractString}
    function Layout{T}(vars::Matrix{T}, titles::Matrix{<:AbstractString}) where T
        size(vars) == size(titles) ||
            throw(ArgumentError("Size of vars and titles must be the same."))
        return new{T}(vars, titles)
    end # function Layout
end # struct Layout
Layout(vars::Matrix{T}, titles::Matrix{<:AbstractString}=string.(vars)) where T = Layout{T}(vars, titles)

Base.size(layout::Layout)::NTuple{2,Int} = size(layout.vars)
Base.axes(layout::Layout, dim::Int)::Base.OneTo{Int64} = axes(layout.vars, dim)
Base.eachindex(layout::Layout)::Base.OneTo{Int64} = eachindex(layout.vars)
Base.getindex(layout::Layout, inx...) = (var=layout.vars[inx...], title=layout.titles[inx...]) # -> @NamedTuple{var::T, title::AbstractString}

struct BackendError <: Exception
    requested::Symbol
end # struct BackendError

function Base.showerror(io::IO, err::BackendError)::Nothing
    if err.requested === :missing
        println(io, "No Makie backend is currently loaded. Please load a backend package first.")
    else # err.requested !== :missing
        println(
            io,
            "Backend package $(err.requested) is not loaded or unsupported. Try `import $(err.requested)` first."
        )
    end # if ===; else
    return nothing
end # function Base.showerror

const classicmodel_layout = Layout(
    [:E  :T  :h],
    [Mk.L"$E$ ($\mathrm{J\,m^{-2}}$)"  Mk.L"$T$ ($\mathrm{\degree\!C}$)"  Mk.L"$h$ ($\mathrm{m}$)"]
)

const mizmodel_layout = Layout(
    [
        :Ew  :Ei  :E
        :Tw  :Ti  :T
        :h   :D   :phi
    ],
    [
        Mk.L"$E_{\mathrm{w}}$ ($\mathrm{J\,m^{-2}}$)"   Mk.L"$E_{\mathrm{i}}$ ($\mathrm{J\,m^{-2}}$)"   Mk.L"$E$ ($\mathrm{J\,m^{-2}}$)"
        Mk.L"$T_{\mathrm{w}}$ ($\mathrm{\degree\!C}$)"  Mk.L"$T_{\mathrm{i}}$ ($\mathrm{\degree\!C}$)"  Mk.L"$T$ ($\mathrm{\degree\!C}$)"
        Mk.L"$\bar{h}$ ($\mathrm{m}$)"                  Mk.L"$\bar{\mathcal{D}}$ ($\mathrm{m}$)"        Mk.L"$\varphi$"
    ]
)

const wimodel_layout = mizmodel_layout

for model in filter(!=(ModelDiff), IU.subtypes(AbstractModel))
    namelower = lowercase(split(string(model), '.')[end])
    @eval default_layout(::$model)::Layout{Symbol} = deepcopy($(Symbol(namelower, "_layout")))
end # for model
function default_layout(::ModelDiff{A,B})::Layout{Symbol} where {A<:AbstractModel, B<:AbstractModel}
    layout = (A === ClassicModel || B === ClassicModel) ?
        deepcopy(classicmodel_layout) : deepcopy(mizmodel_layout)
    foreach(
        i -> layout.titles[i] = Mk.latexstring(raw"$\Delta ", layout.titles[i][2:end]),
        eachindex(layout.titles)
    )
    return layout
end # function default_layout

isloaded(::Val)::Bool = false

init_backend(::Val{S}) where S = throw(BackendError(S))

"""
    backend() -> Union{Module,Missing}

Get the current Makie backend module. If no backend is initialized, returns `missing`.

    backend(bcknd) -> Module

Set the Makie backend to the specified `bcknd` and return the backend module. Supported
backends are `:GLMakie`, `:CairoMakie` and `:WGLMakie`. You need to first load the
corresponding backend package before calling this function.

# Examples
```julia-repl
julia> backend()
missing

julia> import GLMakie; backend(:GLMakie)
GLMakie
```
"""
backend()::Union{Module,Missing} = Mk.current_backend()
backend(bcknd::Symbol)::Module = init_backend(Val(bcknd))

const dlevels = [400(range(0, 400, 20) ./ 400).^2; 450]
# (isD::Val{[Bool]}, diff::Val{[Bool]}, data::Matrix{Float64}) -> (levels, extendlow, extendhigh)
get_levels(::Val{false}, ::Val{false}, _)::Tuple{Int,Nothing,Nothing,Mk.Automatic} = (
    21, nothing, nothing, Mk.automatic
) # normal
get_levels(::Val{true}, ::Val{false}, _)::Tuple{Vector{Float64},Nothing,Symbol,Vector{Int}} = (
    dlevels, nothing, :auto, collect(0:100:400)
) # D sol
get_levels(::Val{true}, ::Val{true}, _)::Tuple{Vector{Float64},Symbol,Symbol,Vector{Int}} = (
    sort!(unique!(d -> isequal(-0.0, d) ? abs(d) : identity(d), [-dlevels; dlevels])),
    :auto, :auto, collect(-400:100:400)
) # D diff
get_levels(::Val{false}, ::Val{true}, data::Matrix{Float64})::Tuple{Vector{Float64},Nothing,Nothing,Mk.Automatic} = (
    maximum(abs, filter(!isnan, data))*range(-1, 1; length=41),
    nothing, nothing, Mk.automatic
) # normal diff

function contourf_tiles(
    t::Vector{T},
    x::Vec,
    layout::Layout{Matrix{Float64}},
    fig::Union{Mk.Figure,Mk.GridPosition};
    xlim::Union{NTuple{2,Real},Nothing}=nothing,
    tlim::Union{NTuple{2,Real}, Nothing}=nothing,
    diff::Bool=false,
    inspect::Bool=false
) where T<:Real # -> Union{Mk.Figure,Mk.GridPosition}
    for row in axes(layout, 1), col in axes(layout, 2)
        subfig = fig[row,col]
        ax = Mk.Axis(
            subfig[1,1];
            title=layout[row,col].title,
            xlabel=(row==lastindex(layout, 1) ? Mk.L"$t$ ($\mathrm{y}$)" : ""),
            ylabel=(col==1 ? Mk.L"x" : ""),
            xticklabelrotation=(T===Int ? 0 : -pi/4),
            limits=(tlim, xlim)
        )
        if all(isnan, layout[row,col].var)
            @warn "All data are NaN at position ($row, $col). Skipping plot."
        else # valid data
            isD = occursin(raw"\mathcal{D}", layout[row,col].title)
            levels, extendlow, extendhigh, ticks = get_levels(Val(isD), Val(diff), layout[row,col].var)
            ctr = Mk.contourf!(ax, t, x, layout[row,col].var; levels, extendlow, extendhigh)
            Mk.Colorbar(subfig[1,2], ctr; ticks)
            if diff
                Mk.contour!(ax, t, x, layout[row,col].var; levels=[0], color=:black)
                cbaxis = Mk.Axis(subfig[1,2]; limits=((-1, 1), (-1, 1)))
                Mk.hidedecorations!(cbaxis)
                Mk.hidespines!(cbaxis)
                Mk.hlines!(cbaxis, 0; color=:black, linewidth=2)
            end # if diff
        end # if all; else
    end # for row, col
    inspect && Mk.DataInspector(fig)
    return fig
end # function contourf_tiles

matricify(vecvec::Vector{Vec})::Matrix{Float64} = permutedims(reduce(hcat, vecvec))

function limit_size(
    xs::Vec, ts::Vector{<:Real},
    xsizelim::Int=1000, tsizelim::Int=1000,
    xrange::NTuple{2,Real}=extrema(xs), trange::NTuple{2,Real}=extrema(ts)
)::NTuple{2,Vector{Int}}
    # find range indices
    tiran = (findfirst(>=(trange[1]), ts), findlast(<=(trange[2]), ts))
    xiran = (findfirst(>=(xrange[1]), xs), findlast(<=(xrange[2]), xs))
    for (iran, name, s) in zip((tiran, xiran), ("time step", "space point"), (ts, xs))
        if any(isnothing, iran) || iran[2]<iran[1] # range ⊄ s
            throw(
                ArgumentError(
                    "No $(name)s stored in the Solutions within the specified range. The range should intersect with $(extrema(s))."
                )
            )
        elseif (iran[1] == iran[2])# only one point in range
            @warn "Only one $name found in the specified range. Nothing will be shown on the contourf plot."
        end # if ||, elseif
    end # for (iran, name, s)
    (xsizelim<=1 || tsizelim<=1) &&
        throw(ArgumentError("Number of points limits must be greater than 1."))
    # limit sizes
    xinx = (xiran[2]-xiran[1]+1) > xsizelim ?
            round.(Int, range(xiran[1], xiran[2], xsizelim)) : # reduce x size
            collect(xiran[1]:xiran[2]) # within the space size limit
    tinx = (tiran[2]-tiran[1]+1) > tsizelim ?
            round.(Int, range(tiran[1], tiran[2], tsizelim)) : # reduce time size
            collect(tiran[1]:tiran[2]) # within the time size limit
    length(tinx)length(xinx) > 1_000_000 &&
        @warn "Number of points to plot $(length(tinx)length(xinx)). This may lead to performance issues."
    return (xinx, tinx)
end # function limit_size

"""
    plot_raw(sols::Solutions, into::Union{Mk.Figure,Mk.GridPosition}=Mk.Figure(); kwargs...) -> Makie.Figure

Plot the the solution variables for each time step in `sols.raw` into the specified Makie
figure or grid position. If `into` is not specified, a new figure will be created. The
function will find available backend if not specified.

# Keyword Arguments
- `layout::Layout{Symbol}`: Layout structure specifying which variables to plot and their
    titles.
- `inspect::Bool`: If true, enables `Makie.DataInspect` for interactive exploration of the
    plot.
- `xsizelim::Int`: Maximum number of spatial points to plot. If the number of spatial
    points in `sols` exceeds this limit, the points will be downsampled uniformly to meet
    the limit.
- `tsizelim::Int`: Maximum number of time steps to plot.
- `xrange::NTuple{2,Real}`: Range of spatial points to plot.
- `trange::NTuple{2,Real}`: Range of time steps to plot.
"""
function plot_raw(
    sols::Solutions{M},
    into::Union{Mk.Figure,Mk.GridPosition}=Mk.Figure();
    layout::Layout{Symbol}=default_layout(M()),
    inspect::Bool=false,
    xsizelim::Int=1000,
    tsizelim::Int=1000,
    xrange::NTuple{2,Real}=extrema(sols.spacetime.x),
    trange::NTuple{2,Real}=extrema(sols.ts)
) where M<:AbstractModel # -> Union{Mk.Figure,Mk.GridPosition}
    xinx, tinx = limit_size(sols.spacetime.x, sols.ts, xsizelim, tsizelim, xrange, trange)
    datatitle = Layout(Matrix{Matrix{Float64}}(undef, size(layout)), layout.titles)
    @simd for linx in eachindex(layout)
        datatitle.vars[linx] = matricify(getindex.(sols.raw[layout[linx].var][tinx], Ref(xinx)))
    end # for inx
    return contourf_tiles(
        sols.ts[tinx], sols.spacetime.x[xinx], datatitle, into;
        xlim=xrange, tlim=trange, inspect, diff=M<:ModelDiff
    )
end # function plot_raw

"""
    plot_avg(sols::Solutions, into::Union{Mk.Figure,Mk.GridPosition}=Mk.Figure(); kwargs...) -> Makie.Figure

Plot the annual average of solution variables in `sols.annual.avg` into the specified Makie
figure or grid position. If `into` is not specified, a new figure will be created. The
function will find available backend if not specified.

# Keyword Arguments
- `layout::Layout{Symbol}`: Layout structure specifying which variables to plot and their
    titles.
- `inspect::Bool`: If true, enables `Makie.DataInspect` for interactive exploration of the
    plot.
- `xsizelim::Int`: Maximum number of spatial points to plot. If the number of spatial
    points in `sols` exceeds this limit, the points will be downsampled uniformly to meet
    the limit.
- `tsizelim::Int`: Maximum number of time steps to plot.
- `xrange::NTuple{2,Real}`: Range of spatial points to plot.
- `trange::NTuple{2,Real}`: Range of time steps to plot.
"""
function plot_avg(
    sols::Solutions{M},
    into::Union{Mk.Figure,Mk.GridPosition}=Mk.Figure();
    layout::Layout{Symbol}=default_layout(M()),
    inspect::Bool=false,
    xsizelim::Int=1000,
    tsizelim::Int=1000,
    xrange::NTuple{2,Real}=extrema(sols.spacetime.x),
    trange::NTuple{2,Real}=(1, sols.spacetime.dur),
) where M<:AbstractModel # -> Union{Mk.Figure,Mk.GridPosition}
    xinx, tinx = limit_size(sols.spacetime.x, collect(1:sols.spacetime.dur), xsizelim, tsizelim, xrange, trange)
    datatitle = Layout(Matrix{Matrix{Float64}}(undef, size(layout)), layout.titles)
    @simd for linx in eachindex(layout)
        datatitle.vars[linx] = matricify(getindex.(sols.annual.avg[layout[linx].var][tinx], Ref(xinx)))
    end # for inx
    return contourf_tiles(
        collect(tinx), sols.spacetime.x[xinx], datatitle, into;
        xlim=xrange, tlim=trange, inspect, diff=M<:ModelDiff
    )
end # function plot_avg

"""
    plot_seasonal(sols::Solutions, fig::Union{Mk.Figure,Mk.GridPosition}=Mk.Figure(); kwargs...) -> Makie.Figure

Using the data from `sols.annual`, plot lines spanned by (`xfunc(sols, year)`,
`yfunc(sols, season, year)`) for each year and for the seasons `:avg`, `:winter`, and
`:summer` into the specified Makie figure or grid position `fig`. If `fig` is not specified,
a new figure will be created. By default, `xfunc` computes the hemispheric mean temperature
from `sols.annual.avg.T`, while `yfunc` computes the ice-covered area using either
concentration `phi` (if it exists) or enthalpy `E`. Lines during the warming period defined
in `sols.forcing` are coloured red, and those during the cooling period are coloured blue.
Lines for the summer peak are dashed, those for winter are thin solid, and those for the
annual average are thick solid.

# Keyword Arguments
- `xfunc::Function`: A function that takes in `sols` and `year` and returns a `Float64`
    representing the x-coordinate for that year.
- `yfunc::Function`: A function that takes in `sols`, `season`, and `year` and returns a
    `Float64` representing the y-coordinate for that season and year.
- `title::AbstractString`: Title of the plot.
- `xlabel::AbstractString`: Label for the x-axis.
- `ylabel::AbstractString`: Label for the y-axis.
- `inspect::Bool`: If true, enables `Makie.DataInspect` for interactive exploration of the
    plot.
"""
function plot_seasonal(
    sols::Solutions{<:AbstractModel,F,true},
    fig::Union{Mk.Figure,Mk.GridPosition}=Mk.Figure();
    xfunc::Function=((sols, year) -> hemispheric_mean(sols.annual.avg.T[year], sols.spacetime.x)),
    yfunc::Function=ice_area,
    title::AbstractString="Ice covered area",
    xlabel::AbstractString=Mk.L"$\tilde{T}$ ($\mathrm{\degree\!C}$)",
    ylabel::AbstractString=Mk.L"A_i",
    inspect::Bool=false
) where F # -> Union{Mk.Figure,Mk.GridPosition}
    xdata = xfunc.(Ref(sols), 1:sols.spacetime.dur)
    ax = Mk.Axis(fig[1,1]; title, xlabel, ylabel)
    groups = (
        Warming=Vector{Mk.Lines{Tuple{Vector{Mk.Point{2,Float64}}}}}(),
        Cooling=Vector{Mk.Lines{Tuple{Vector{Mk.Point{2,Float64}}}}}()
    )
    for (domain, group, inx, colour) in zip(
        keys(groups),
        values(groups),
        (sols.forcing.domain[2]:sols.forcing.domain[3], sols.forcing.domain[4]:sols.forcing.domain[5]),
        (Mk.Cycled(6), Mk.Cycled(1))
    ), season in (:avg, :winter, :summer)
        width = 1
        season === :avg && (width += (domain===:Warming ? 2 : 1))
        push!(
            group,
            Mk.lines!(
                ax, xdata[inx], yfunc.(Ref(sols), Ref(season), inx);
                color=colour, linewidth=width, linestyle=(season===:summer ? :dash : :solid)
            )
        ) # push!
    end # for domain, inx, colour, season
    Mk.Legend(
        fig[1,2],
        collect(values(groups)),
        fill(["mean", "winter", "summer"], 2),
        string.(collect(keys(groups)))
    )
    inspect && Mk.DataInspector(fig)
    return fig
end # function plot_seasonal

import PrecompileTools as PT

function precompile(bcknd::Module)::Nothing
    PT.@setup_workload begin
        ints = collect(1:10)
        floats = collect(0.1:0.1:1.0)
        x = collect(0.1:0.1:1.0)
        layout = Layout(
            reshape([rand(10, 10)], 1, 1), reshape(AbstractString[Mk.L"title"], 1, 1)
        )
        bcknd.activate!()
        PT.@compile_workload begin
            for t in (ints, floats)
                contourf_tiles(t, x, layout, Mk.Figure(); inspect=true)
            end # for t
        end # PT.@compile_workload begin
    end # PT.@setup_workload begin
    return nothing
end # function precompile

end # module Plot
