module Plot

using ..Utilities, ..Infrastructure

import Makie

export Layout, backend
export plot_raw, plot_avg, plot_seasonal

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

function init_backend(val::Val)
    name = typeof(val).parameters[1]
    if val isa Val{:GLMakie} || val isa Val{:CairoMakie}
        throw(ArgumentError("Backend not loaded. Please load the backend package $name first."))
    else
        throw(ArgumentError("Unsupported backend $name."))
    end # if ||, else
end

backend()::Union{Module,Missing} = Makie.current_backend()
backend(bcknd::Symbol)::Module = init_backend(Val(bcknd))

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
        ctr = Makie.contourf!(ax, t, x, layout[row,col].var)
        Makie.Colorbar(subfig[1,2], ctr)
    end # for row, col
    @eval Main innerfig = $fig
    return fig
end # function contourf_tiles

matricify(vecvec::Vector{Vec})::Matrix{Float64} = permutedims(reduce(hcat, vecvec))

function plot_raw(
    sols::Solutions{F,C},
    bcknd::Symbol=:GLMakie;
    layout::Layout{Symbol}=(:phi in propertynames(sols.raw) ? miz_layout : classic_layout)
)::Makie.Figure where {F, C}
    backend(bcknd)
    datatitle = Layout(Matrix{Matrix{Float64}}(undef, size(layout)), layout.titles)
    @simd for inx in eachindex(layout)
        datatitle.vars[inx] = matricify(getproperty(sols.raw, layout[inx].var))
    end # for inx
    return contourf_tiles(sols.ts, sols.spacetime.x, datatitle)
end # function plot_raw

function plot_avg(
    sols::Solutions{F,C},
    bcknd::Symbol=:GLMakie;
    layout::Layout{Symbol}=(:phi in propertynames(sols.raw) ? miz_layout : classic_layout)
)::Makie.Figure where {F, C}
    backend(bcknd)
    datatitle = Layout(Matrix{Matrix{Float64}}(undef, size(layout)), layout.titles)
    @simd for inx in eachindex(layout)
        datatitle.vars[inx] = matricify(getproperty(sols.seasonal.avg, layout[inx].var))
    end # for inx
    return contourf_tiles(collect(1:sols.spacetime.dur), sols.spacetime.x, datatitle)
end # function plot_avg

function plot_seasonal(
    sols::Solutions{F,false},
    bcknd::Symbol=:GLMakie;
    xfunc::Function=(
        (sols::Solutions{F,false}, year::Int) ->
            hemispheric_mean(sols.seasonal.avg.T[year], sols.spacetime.x)
    ),
    yfunc::Function=(
        :phi in propertynames(sols.raw) ?
            (
                (sols::Solutions{F,false}, season::Symbol, year::Int) ->
                    2.0*pi * hemispheric_mean(getproperty(sols.seasonal, season).phi[year], sols.spacetime.x)
            ) :
            (
                (sols, season, year) ->
                    2.0*pi * hemispheric_mean((getproperty(sols.seasonal, season).E[year]<0.0), sols.spacetime.x)
            ) # ? :
    ), # ()
    title::AbstractString="Ice covered area",
    xlabel::AbstractString=Makie.L"$\tilde{\mathsf{T}}$ ($\mathrm{\degree\!C}$)",
    ylabel::AbstractString=Makie.L"A_i$"
)::Makie.Figure where F
    backend(bcknd)
    xdata = xfunc.(Ref(sols), sols.spacetime.dur)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1,1]; title=title, xlabel=xlabel, ylabel=ylabel)
    groups = (
        Warming=Vector{Makie.Lines{Tuple{Vector{Makie.Point{2,Float64}}}}}(),
        Cooling=Vector{Makie.Lines{Tuple{Vector{Makie.Point{2,Float64}}}}}()
    )
    for (domain, group, inx, colour) in zip(
        keys(groups),
        values(groups),
        (sols.forcing.domain[2]:sols.forcing.domain[3], sols.forcing.domain[4]:sols.forcing.domain[5]),
        (Makie.wong_colors()[6], Makie.wong_colors()[1])
    ), season in (:avg, :winter, :summer)
        width = 1.0
        if season === :avg
            width += domain===:Warming ? 2.0 : 1.0
        end # if ===
        push!(
            group,
            Makie.lines!(
                ax, xdata[inx], yfunc.(Ref(sols), Ref(season), sols.spacetime.dur[inx]);
                color=colour, linewidth=width, linestyle=(season===:summer ? :dash : :solid)
            )
        ) # push!
    end # for domain, inx, colour, season
    Makie.Legend(
        fig[1,2], values(groups), fill(["mean", "winter", "summer"], 2), keys(groups)
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
