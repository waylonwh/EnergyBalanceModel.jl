using EnergyBalanceModel, Test

# shared helpers, fixtures and lazily integrated solutions
include("testsetup.jl")

@testset verbose = true "EnergyBalanceModel.jl" begin
    include("test_utilities.jl")
    include("test_infrastructure.jl")
    include("test_classicebm.jl")
    include("test_mizebm.jl")
    include("test_wimebm.jl")
    include("test_solvers.jl")
    include("test_plot.jl")
    include("test_integration.jl")
end # @testset "EnergyBalanceModel.jl"
