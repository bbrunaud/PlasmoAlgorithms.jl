using JuMP
using CPLEX
using Plasmo
using PlasmoAlgorithms
include("input.jl")
include("master.jl")
include("subproblem.jl")
#define master problem in financeplanning 
master = generate_master()

g = ModelGraph()
setsolver(g, CplexSolver())
n1 = add_node(g)
setmodel(g, master)










