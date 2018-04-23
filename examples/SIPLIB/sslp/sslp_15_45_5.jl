using Gurobi
include("../generateSIPLIB.jl")

const PA = PlasmoAlgorithms
g = genproblem("sslp",15,45,5)
g.solver = GurobiSolver(MIPGap=0.01)
