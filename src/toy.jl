using Plasmo 
using JuMP
using CPLEX
import Plasmo.solve, Base.==
using Requires
using PyCall
using CPLEX
using MathProgBase
master = Model(solver=CplexSolver())
@variable(master, 0<=x1<=1)
@variable(master, x2, Bin)

@constraint(master, -x1 -x2 >= -1.5)
@constraint(master, l1, x1 == 0.3)
@constraint(master, l2, x2==0)
sub = Model(solver=CplexSolver()) 
@variable(sub, 0<=y1<=1)
@variable(sub, 0<=y2<=1)
@variable(sub, y3, Bin)
@variable(sub, y4, Bin)
@variable(sub, 0<=x1<=1)
@variable(sub, x2, Bin)

@constraint(sub, -2*y1-3*y2-4*y3-5*y4>= -5-0.3*x1)
@constraint(sub, -6*y1-y2-3*y3-2*y4>=-10+0.3*x2)
@objective(sub, Min, -16*y1-19*y2-23*y3-28*y4)
sub1 = deepcopy(sub)
@objective(sub1, Min, (x1-0.5)^2 + (x2-0.5)^2)
g = PlasmoGraph()
g.solver = CplexSolver()
n1 = add_node(g)
setmodel(n1, master)
n2 = add_node(g)
setmodel(n2, sub)
@linkconstraint(g, n1[:x1]==n2[:x1])
@linkconstraint(g, n1[:x2]==n2[:x2])
edge = add_edge(g, n1, n2)

