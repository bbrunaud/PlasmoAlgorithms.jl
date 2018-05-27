using PlasmoAlgorithms
using Plasmo
using Gurobi

g = smpsread("dcap/dcap233_200")
g.solver = GurobiSolver(OutputFlag=0)
temp_g = deepcopy(g)
bendersolve(g, cuts=[:Root], max_iterations=1)
g2 = PlasmoGraph()
n1 = add_node(g2,getnodes(g)[2].model)
n2 = add_node(g2,getnodes(g)[1].model)
# n3 = add_node(g2,getnodes(g)[3].model)
edge1 = add_edge(g2, n1, n2)
# edge2 = add_edge(g2, n1, n3)
y = getindex(n1.model,:y)
@linkconstraint(g2, [i in keys(y)], n1[:y][i] == n2[:x][i])
# @linkconstraint(g2, [i in keys(y)], n1[:y][i] == n3[:x][i])
g2.solver = GurobiSolver()
# bendersolve(g2, cuts=[:Root], max_iterations=1)