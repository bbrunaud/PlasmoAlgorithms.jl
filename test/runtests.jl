using Base.Test

@testset "Lagrange" begin include("test_lagrange.jl") end
@testset "Benders" begin include("test_benders.jl") end
