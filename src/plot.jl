module Plot

using ..Utilities, ..Infrastructure

import Makie

export Layout, backend
export plot_raw, plot_avg, plot_seasonal

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
    function Layout{T}(vars::Matrix{T}, titles::Matrix{AbstractString}) where T
        if size(vars) != size(titles)
            throw(ArgumentError("Size of vars and titles must be the same."))
        end # if !=
        return new{T}(vars, titles)
    end # function Layout
    Layout(vars::Matrix{T}, titles::Matrix{AbstractString}) where T = Layout{T}(vars, titles)
end

(Base.size(layout::Layout{T})::NTuple{2,Int}) where T = size(layout.vars)
(Base.axes(layout::Layout{T}, dim::Int)::Base.OneTo{Int64}) where T = axes(layout.vars, dim)
(Base.eachindex(layout::Layout{T})::Base.OneTo{Int64}) where T = eachindex(layout.vars)
(Base.getindex(layout::Layout{T}, inx...)::@NamedTuple{var::T, title::AbstractString}) where T =
    (var=layout.vars[inx...], title=layout.titles[inx...])

const miz_layout = Layout(
    [
        :Ew :Ei :E
        :Tw :Ti :T
        :h  :D  :phi
    ],
    AbstractString[
        Makie.L"$E_w$ ($\mathrm{J\,m^{-2}}$)"  Makie.L"$E_i$ ($\mathrm{J\,m^{-2}}$)"       Makie.L"$E$ ($\mathrm{J\,m^{-2}}$)"
        Makie.L"$T_w$ ($\mathrm{\degree\!C}$)" Makie.L"$T_i$ ($\mathrm{\degree\!C}$)"      Makie.L"$T$ ($\mathrm{\degree\!C}$)"
        Makie.L"$\bar{h}$ ($\mathrm{m}$)"      Makie.L"$\bar{\mathcal{D}}$ ($\mathrm{m}$)" Makie.L"\varphi"
    ]
)

const classic_layout = Layout(
    [:E :T :h],
    AbstractString[Makie.L"$E$ ($\mathrm{J\,m^{-2}}$)" Makie.L"$T$ ($\mathrm{\degree\!C}$)" Makie.L"$h$ ($\mathrm{m}$)"]
)

default_layout(::MIZModel) = miz_layout
default_layout(::ClassicModel) = classic_layout

isloaded(::Val)::Bool = false

function find_backend()::Union{Symbol,Nothing}
    for backend in (:GLMakie, :CairoMakie, :WGLMakie)
        if isloaded(Val(backend))
            return backend
        end # if isloaded
    end # for backend
    return nothing
end # function find_backend

function init_backend(val::Val) # -> ERROR
    name = typeof(val).parameters[1]
    loaded = find_backend()
    errmsg = "Backend package $name is not loaded or unsupported. Please load the the backend package first."
    if !isnothing(loaded)
        errmsg *= "\nHint: Another backend package $loaded is already loaded."
    end # if !
    throw(ArgumentError(errmsg))
end # function init_backend

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
backend()::Union{Module,Missing} = Makie.current_backend()
backend(bcknd)::Module = init_backend(Val(bcknd))

function contourf_tiles(t::Vector{T}, x::Vec, layout::Layout{Matrix{Float64}})::Makie.Figure where T<:Real
    fig = Makie.Figure()
    for row in axes(layout, 1), col in axes(layout, 2)
        subfig = fig[row,col]
        ax = Makie.Axis(
            subfig[1,1];
            title=layout[row,col].title,
            xlabel=(row==lastindex(layout, 1) ? Makie.L"$t$ ($\mathrm{y}$)" : ""),
            ylabel=(col==1 ? Makie.L"x" : ""),
            limits=(nothing, (0, 1))
        )
        if all(isnan, layout[row,col].var)
            @warn "All data are NaN at position ($row, $col). Skipping plot."
        else # valid data
            ctr = Makie.contourf!(ax, t, x, layout[row,col].var)
            Makie.Colorbar(subfig[1,2], ctr)
        end # if all; else
    end # for row, col
    return fig
end # function contourf_tiles

matricify(vecvec::Vector{Vec})::Matrix{Float64} = permutedims(reduce(hcat, vecvec))

"""
    plot_raw(sols::Solutions{M,F,C},bcknd::Union{Symbol,Nothing}=...; layout::Layout{Symbol}=... -> Makie.Figure

Plot the the solution variables for each time step in `sols.raw` using the specified Makie
backend `bcknd` and `layout`. The function will find available backend if not specified. By
default, the layout is set to `miz_layout` if sols is a `Solutions{MIZModel}`, and
`classic_layout` if sols is a `Solutions{ClassicModel}`. Use
`EnergyBalanceModel.Plot.default_layout(miz)` or
`EnergyBalanceModel.Plot.default_layout(classic)` to get default layouts.

"""
function plot_raw(
    sols::Solutions{M,F,C},
    bcknd::Union{Symbol,Nothing}=find_backend();
    layout::Layout{Symbol}=default_layout(M())
)::Makie.Figure where {M<:AbstractModel, F, C}
    backend(bcknd)
    datatitle = Layout(Matrix{Matrix{Float64}}(undef, size(layout)), layout.titles)
    @simd for inx in eachindex(layout)
        datatitle.vars[inx] = matricify(getproperty(sols.raw, layout[inx].var))
    end # for inx
    return contourf_tiles(sols.ts, sols.spacetime.x, datatitle)
end # function plot_raw

"""
    plot_avg(sols::Solutions{<:AbstractModel,F,C}, bcknd::Union{Symbol,Nothing}=...; layout::Layout{Symbol}=... -> Makie.Figure

Plot the annual average of solution variables in `sols.annual.avg` using the specified
Makie backend `bcknd` and `layout`.
"""
function plot_avg(
    sols::Solutions{M,F,C},
    bcknd::Union{Symbol,Nothing}=find_backend();
    layout::Layout{Symbol}=default_layout(M())
)::Makie.Figure where {M<:AbstractModel, F, C}
    backend(bcknd)
    datatitle = Layout(Matrix{Matrix{Float64}}(undef, size(layout)), layout.titles)
    @simd for inx in eachindex(layout)
        datatitle.vars[inx] = matricify(getproperty(sols.annual.avg, layout[inx].var))
    end # for inx
    return contourf_tiles(collect(1:sols.spacetime.dur), sols.spacetime.x, datatitle)
end # function plot_avg

(ice_area(sols::Solutions{ClassicModel,F,C}, season::Symbol, year::Int)::Float64) where {F, C} =
    2.0pi * hemispheric_mean((getproperty(sols.annual, season).E[year].<0.0), sols.spacetime.x)
(ice_area(sols::Solutions{MIZModel,F,C}, season::Symbol, year::Int)::Float64) where {F, C} =
    2.0pi * hemispheric_mean(getproperty(sols.annual, season).phi[year], sols.spacetime.x)

"""
    plot_seasonal(sols::Solutions{F,false}, bcknd::Union{Symbol,Nothing}=...; kwargs...) -> Makie.Figure

Using the data from `sols.annual`, plot lines spanned by (`xfunc(sols, year)`,
`yfunc(sols, season, year)`) for each year and for the seasons `:avg`, `:winter`, and
`:summer`. By default, `xfunc` computes the hemispheric mean temperature from
`sols.annual.avg.T`, while `yfunc` computes the ice-covered area using either
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
"""
function plot_seasonal(
    sols::Solutions{M,F,false},
    bcknd::Union{Symbol,Nothing}=find_backend();
    xfunc::Function=((sols, year) -> hemispheric_mean(sols.annual.avg.T[year], sols.spacetime.x)),
    yfunc::Function=ice_area,
    title::AbstractString="Ice covered area",
    xlabel::AbstractString=Makie.L"$\tilde{T}$ ($\mathrm{\degree\!C}$)",
    ylabel::AbstractString=Makie.L"A_i"
)::Makie.Figure where {M<:AbstractModel, F}
    backend(bcknd)
    xdata = xfunc.(Ref(sols), 1:sols.spacetime.dur)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1,1]; title, xlabel, ylabel)
    groups = (
        Warming=Vector{Makie.Lines{Tuple{Vector{Makie.Point{2,Float64}}}}}(),
        Cooling=Vector{Makie.Lines{Tuple{Vector{Makie.Point{2,Float64}}}}}()
    )
    for (domain, group, inx, colour) in zip(
        keys(groups),
        values(groups),
        (sols.forcing.domain[2]:sols.forcing.domain[3], sols.forcing.domain[4]:sols.forcing.domain[5]),
        (Makie.Cycled(6), Makie.Cycled(1))
    ), season in (:avg, :winter, :summer)
        width = 1.0
        if season === :avg
            width += domain===:Warming ? 2.0 : 1.0
        end # if ===
        push!(
            group,
            Makie.lines!(
                ax, xdata[inx], yfunc.(Ref(sols), Ref(season), inx);
                color=colour, linewidth=width, linestyle=(season===:summer ? :dash : :solid)
            )
        ) # push!
    end # for domain, inx, colour, season
    Makie.Legend(
        fig[1,2],
        collect(values(groups)),
        fill(["mean", "winter", "summer"], 2),
        string.(collect(keys(groups)))
    )
    return fig
end # function plot_seasonal

function unsafesave(plt::Makie.Figure, path::String; spwarn::Bool=false, kwargs...)::String
    if !spwarn
        @warn "`unsafesave` may overwrite existing files. Use `save` instead."
    end # if !
    Makie.save(path, plt; kwargs...)
    return path
end

end # module Plot
