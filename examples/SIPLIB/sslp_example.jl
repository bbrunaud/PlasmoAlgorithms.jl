
using PlasmoAlgorithms
using Plasmo
using Gurobi

g = smpsread("sslp/sslp_15_45_5")
n = collect(values(getnodes(g)))

g.solver = GurobiSolver()

# n[2] is the root node