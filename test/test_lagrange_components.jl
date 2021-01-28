using PlasmoAlgorithms
using Test

const PA = PlasmoAlgorithms

include("fisher.jl")

@test typeof(graph) == OptiGraph

node1 = getnode(graph, 1)
@test typeof(node1) == OptiNode

# Prepare Function
@test PA.lgprepare(graph,  δ=0.5, maxnoimprove=2,cpbound=1e6)
m1 = getmodel(node1)
@test objective_sense(m1) == MOI.MIN_SENSE
@test PA.getattribute(graph, :numlinks) == 2
@test PA.getattribute(graph, :normalized) == -1

# Inital Relaxation
λ0 = PA.getattribute(graph,:λ)[end]
@test λ0 == zeros(2)
lpobj = PA.initialrelaxation(graph)
@test lpobj == -18.0
λlp = PA.getattribute(graph,:λ)[end]
@test λlp == [8.0, 2.0]

nodes = [node for node in collect(getnodes(graph))]
@test length(nodes) == 2
