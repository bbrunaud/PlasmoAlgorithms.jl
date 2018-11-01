using JuMP
using CPLEX
using Plasmo
using PlasmoAlgorithms

mp = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0))
sp1 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0))
sp2 = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0))


@variable(mp,x1, Bin)
@objective(mp,Min,-x1)

@variable(sp1, y1, Bin)
@variable(sp1, y2>=0)
@variable(sp1, x1, Bin)
@constraint(sp1, -3*y1+y2>=-6 + 4*x1)
@objective(sp1,Min, -2*y1 + 4 * y2)

@variable(sp2, y1, Bin)
@variable(sp2, y2>=0)
@variable(sp2, x1, Bin)
@constraint(sp2, -3*y1+y2>=-2 + 4*x1)
@objective(sp2, Min, -2*y1 + 4 * y2)



g = ModelGraph()
setsolver(g, CplexSolver())
n1 = add_node(g)
n2 = add_node(g)
n3 = add_node(g)


setmodel(n1, mp)
setmodel(n2, sp1)
setmodel(n3, sp2)


edge1 = Plasmo.add_edge(g,n1,n2)
edge2 = Plasmo.add_edge(g,n1,n3)



@linkconstraint(g, n1[:x1] == n2[:x1])
@linkconstraint(g, n1[:x1] == n3[:x1])


PlasmoAlgorithms.crosssolve(g,max_iterations = 10, bdcuts =[:LP])
