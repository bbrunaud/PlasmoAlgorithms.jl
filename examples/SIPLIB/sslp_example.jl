using PlasmoAlgorithms
using Plasmo
using CPLEX

g = smpsread("sslp/sslp_15_45_15")
g.solver = CplexSolver(CPX_PARAM_SCRIND=0)
# g1 = deepcopy(g)
# bendersolve(g, cuts=[:LP], max_iterations=100, timelimit=10000)
bendersolve(g, cuts=[:GMI], max_iterations=60, timelimit=100000)
# global optimal_solutions = Dict()
# all_nodes= getnodes(g)

# for i in 1:length(all_nodes)
# 	optimal_solutions[all_nodes[i].label] = all_nodes[i].model.colVal
# end
# g1.attributes[:optimal_solutions] = optimal_solutions
# bendersolve(g1, cuts=[:Root], max_iterations=20)
