module EnergyBalanceModel

import SparseArrays

Parameters = Dict{Symbol, Float64}
ParamSet = Set{Symbol}

struct SpaceTime
    nx::Int # number of evenly spaced latitudinal gridboxes (equator to pole)
    dx::Float64 # grid box width
    x::Vector{Float64} # grid
    nt::Int # number of timesteps per year (limited by numerical stability)
    dt::Float64 # timestep
    t::Vector{Float64} # time vector in a year

    function SpaceTime(nx::Int, nt::Int)
        dx = 1.0 / nx
        x = collect(range(dx/2.0, 1.0 - dx/2.0, nx))
        dt = 1.0 / nt
        t = collect(range(dt/2.0, 1.0 - dt/2.0, nt))
        return new(nx, dx, x, nt, dt, t)
    end
end

struct Solution
    dur::Int # duration of simulation in years
    lastonly::Bool # store only last year of solution
    T::Vector{Float64} # time vector for solution storage
    sol::Dict{Symbol, Vector{Vector{Float64}}} # solution storage

    function Solution(
        duration::Int, lastonly::Bool, st::SpaceTime, vars::Set{Symbol}; debug::Bool=false
    )
        dur_store = lastonly ? 1 : duration # for calculating the number of total timesteps
        T = st.dt : st.dt : dur_store
        debug ? push!(vars, :debug) : nothing
        foreach(var -> (sol[var] = Vector{Vector{Float64}}(undef, length(T))), vars)
    end
end

default_params = Parameters(
    :D => 0.6, # diffusivity for heat transport (W m^-2 K^-1)
    :A => 193.0, # OLR when T = T_m (W m^-2)
    :B => 2.1, # OLR temperature dependence (W m^-2 K^-1)
    :cw => 9.8, # ocean mixed layer heat capacity (W yr m^-2 K^-1)
    :S0 => 420.0, # insolation at equator  (W m^-2)
    :S1 => 338.0, # insolation seasonal dependence (W m^-2)
    :S2 => 240.0, # insolation spatial dependence (W m^-2)
    :a0 => 0.7, # ice-free co-albedo at equator
    :a2 => 0.1, # ice-free co-albedo spatial dependence
    :ai => 0.4, # co-albedo where there is sea ice
    :Fb => 4.0, # heat flux from ocean below (W m^-2)
    :k => 2.0, # sea ice thermal conductivity (W m^-2 K^-1)
    :Lf => 9.5, # sea ice latent heat of fusion (W yr m^-3)
    :F => 0.0, # radiative forcing (W m^-2)
    :cg => 0.01 * 9.8, # ghost layer heat capacity(W yr m^-2 K^-1)
    :tau => 1e-5, # ghost layer coupling timescale (yr)
    :Tm => 0.0, # mean temperature (C)
    :m1 => 1.6e-6 * 31536000, # empirical constants of lateral melt
    :m2 => 1.36, # empirical constants of lateral melt
    :alpha => 0.66, # floe geometry constant, Ai = alpha * D^2
    :rl => 0.5, # lead region width (m)
    :Dmin => 1.0, # new pancake size (m)
    :Dmax => 156, # largest floe length (m)
    :hmin => 0.1, # new pancake thickness (m)
    :kappa => 0.01  * 31536000 # floe welding parameter
)

# Create a parameter dictionary from default values for a given Set
function get_defaultpar(paramset::ParamSet)::Parameters
    setvec = collect(paramset)
    return Parameters(setvec .=> getindex.(Ref(default_params), setvec))
end

# calculate diffusion operator matrix
function get_diffop!(par::Parameters, st::SpaceTime)::Nothing
    xb = st.dx : st.dx : 1.0 - st.dx
    lambda = @. par[:D]/st.dx^2 * (1 - xb^2)
    l1 = pushfirst!(-copy(lambda), 0.0)
    l2 = push!(-copy(lambda), 0.0)
    l3 = -l1 - l2
    par[:diffop] = SparseArrays.spdiagm(-1 => -l1[2:st.nx], 0 => -l3, 1 => -l2[1:st.nx-1]) # !
    return nothing
end

diffusion(f::Vector{Float64}, diffop::AbstractMatrix{Float64})::Vector{Float64} = diffop * f
D∆ = diffusion

module MIZEBM

import ..Parameters, ..ParamSet

miz_paramset = ParamSet(
    [
        :D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :Tm, :m1, :m2, :alpha,
        :rl, :Dmin, :Dmax, :hmin, :kappa
    ]
)

# Functions of space and time
solar(x::Vector{Float64}, t::Float64, par::Parameters)::Vector{Float64} =
    @. par[:S0] - par[:S1] * x * cos(2 * pi * t) - par[:S2] * x^2
coalbedo(x::Vector{Float64}, ice::Bool, par::Parameters)::Vector{Float64} =
    ice ? fill(par[:ai], length(x)) : @. par[:a0] - par[:a2] * x^2

# temperatures
water_temp(Ew::Vector{Float64}, phi::Vector{Float64}, par::Parameters)::Vector{Float64} =
    @. par[:Tm] + Ew / ((1-phi) * par[:cw])
function ice_temp(
    x::Vector{Float64},
    t::Float64,
    h::Vector{Float64},
    Tw::Vector{Float64},
    phi::Vector{Float64},
    f::Float64,
    par::Parameters
)::Vector{Float64}
    return nothing # placeholder
end # module MIZEBM
function Tbar(Ti::Vector{Float64}, Tw::Vector{Float64}, phi::Vector{Float64})::Vector{Float64}
    return nothing # placeholder
end


module ClassicEBM

import ..Parameters, ..ParamSet

classic_paramset = ParamSet(
    [:D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :F, :cg, :tau]
)

end # module ClassicEBM


end # module EnergyBalanceModel
