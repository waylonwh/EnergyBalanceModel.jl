module DistributedExt

using EnergyBalanceModel
import EnergyBalanceModel as EBM

import EnergyBalanceModel.Infrastructure: Model

import Distributed as Distb

struct Problem{M<:Model,F,C}
    model::M
    spacetime::SpaceTime{F}
    forcing::Forcing{C}
    parameters::Collection{Float64}
    initconds::Collection{Vec}
end

function cartesian_problems(
    models::Vector{Model}, sts::Vector{SpaceTime}, forcings::Vector{Forcing},
    pars::Vector{Collection{Float64}}, inits::Vector{Collection{Vec}};
    lastonlies::Vector{Bool}=[true]
)::Vector{Problem}
    problems = Vector{Problem}()
    for pt in Iterators.product(models, sts, forcings, pars, inits, lastonlies)
        push!(problems, Problem(pt...))
    end # for pt
    return problems
end # function cartesian_problems

function needworkers(n::Int)::Vector{Int}
    if Distb.nprocs() - 1 < n # add workers if not enough
        added = Distb.addprocs(n - Distb.nprocs() + 1)
        Distb.@everywhere added using EnergyBalanceModel
    end # if <
    return Distb.workers()
end # function needworkers

# Add methods to EBM.Infrastructure.integrate
(EBM.Infrastructure.integrate(problem::Problem{M,F,C}, kwargs...)::Solutions{M,F,C}) where {M<:Model, F, C} = integrate(
    problem.model, problem.spacetime, problem.forcing, problem.parameters, problem.initconds;
    kwargs...
)

function EBM.Infrastructure.integrate(
    problems::Vector{Problem}, runon::Vector{Int}=Distb.workers();
    progress::Bool=true, kwargs...
)::Vector{Solutions}

end # function EBM.Infrastructure.integrate

# TODO
function EBM.Infrastructure.integrate(
    problems::Vector{Problem}, nworkers::Int; kwargs...
)::Vector{Solutions} end # function EBM.Infrastructure.integrate

end # module DistributedExt


import Distributed as Distb
