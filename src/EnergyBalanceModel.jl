module EnergyBalanceModel

import SparseArrays

# The set of parameters applicable to a given model
struct ParameterSet <: AbstractSet
    set::Set{Symbol}
end

struct SpaceTime
    nx::Int # number of evenly spaced latitudinal gridboxes (equator to pole)
    dx::Float64 # grid box width
    x::Vector{Float64} # grid
    nt::Int # number of timesteps per year (limited by numerical stability)
    dt::Float64 # timestep
    t::Vector{Float64} # time vector
    function SpaceTime(nx::Int, nt::Int)
        dx = 1.0 / nx
        x = collect(range(dx/2.0, 1.0 - dx/2.0, nx))
        dt = 1.0 / nt
        t = collect(range(dt/2.0, 1.0 - dt/2.0, nt))
        return new(nx, dx, x, nt, dt, t)
    end
end

# Structure to store parameters
struct Parameters <: AbstractDict{Symbol, Float64}
    paramval::Dict{Symbol, Float64}
end

Base.getindex(p::Parameters, key::Symbol) = p.paramval[key]
Base.setindex!(p::Parameters, val::Float64, key::Symbol) = setindex!(p.paramval, val, key)

classic_paramset = ParameterSet(
    Set([:D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :F, :cg, :tau])
)
miz_paramset = ParameterSet(
    Set(
        [
            :D, :A, :B, :cw, :S0, :S1, :S2, :a0, :a2, :ai, :Fb, :k, :Lf, :Tm, :m1, :m2,
            :alpha, :rl, :Dmin, :Dmax, :hmin, :kappa
        ]
    )
)

all_paramset = ParameterSet(union(classic_paramset.set, miz_paramset.set))
default_params = Parameters(
    Dict(
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
)

# Create a Parameters object from default values for a given ParameterSet
function Parameters(paramset::ParameterSet)
    setvec = collect(paramset.set)
    return Parameters(Dict(setvec .=> getindex.(Ref(default_params.paramval), setvec)))
end

# calculate diffusion operator matrix
function get_diffop!(par::Parameters, st::SpaceTime)::Nothing
    xb = st.dx : st.dx : 1.0 - st.dx
    lambda = par[:D]/st.dx^2 * (1 .- xb.^2)
    l1 = pushfirst!(-copy(lambda), 0.0)
    l2 = push!(-copy(lambda), 0.0)
    l3 = -l1 - l2
    par[:diffop] = SparseArrays.spdiagm(-1 => -l1[2:st.nx], 0 => -l3, 1 => -l2[1:st.nx-1])
    return nothing
end

module ClassicEBM

end # module ClassicEBM


end # module EnergyBalanceModel
