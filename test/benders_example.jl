using JuMP
using Xpress
using Plasmo
using PlasmoAlgorithms

graph = OptiGraph()
optimizer  = Xpress.Optimizer

##Place MP and SP into PlasmoGraph
mp = @optinode(graph)
JuMP.set_optimizer(getmodel(mp), optimizer)
@variable(mp,y>=0)
@objective(mp,Min,2y)

sp = @optinode(graph)
JuMP.set_optimizer(getmodel(sp), optimizer)
@variable(sp,x[1:2]>=0)
@variable(sp,y>=0)
@constraint(sp,c1,2x[1]-x[2]+3y>=4)
@constraint(sp,x[1]+2x[2]+y>=3)
@objective(sp,Min,2x[1]+3x[2])


##Set sp as a child node of ms
edge = add_edge!(graph,[mp,sp])

## Linking constraints between MP and SP
@linkconstraint(graph, mp[:y] == sp[:y])
Plasmo.set_optimizer(graph, optimizer)
