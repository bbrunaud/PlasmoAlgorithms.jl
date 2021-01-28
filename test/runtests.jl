using Test

@testset "Lagrange components" begin include("test_lagrange_components.jl") end
@testset "Lagrange Subgradient" begin include("test_lagrange_subgradient.jl") end
@testset "Lagrange Cutting Planes" begin include("test_lagrange_cutting_planes.jl") end
@testset "Benders" begin include("test_benders.jl") end