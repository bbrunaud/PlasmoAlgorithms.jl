using PlasmoAlgorithms
using Plasmo
using Gurobi

g = smpsread("dcap/dcap233_200")
n = collect(values(getnodes(g)))

g.solver = GurobiSolver()