using PlasmoAlgorithms
using Plasmo
using Gurobi

g = smpsread("sslp/sslp_5_25_50")
g.solver = GurobiSolver(OutputFlag=0)
g1 = deepcopy(g)
bendersolve(g, cuts=[:LP], max_iterations=100)
global optimal_solutions = Dict()
all_nodes= getnodes(g)

for i in 1:length(all_nodes)
	optimal_solutions[all_nodes[i].label] = all_nodes[i].model.colVal
end
g1.attributes[:optimal_solutions] = optimal_solutions
bendersolve(g1, cuts=[:Root], max_iterations=20)
