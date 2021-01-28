using Test
using PlasmoAlgorithms

const PA = PlasmoAlgorithms

include("benders_example.jl")

@test PA.identifylevels(graph) == 2

@test PA.bdprepare(graph)

bd = bendersoptimize!(graph, verbose=false)

@test bd.objval ≈ 5.285 atol=0.01
@test bd.bestbound ≈ 5.285 atol=0.01
@test bd.gap == 0


