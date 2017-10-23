using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi
using DataFrames

Logging.configure(level=DEBUG)

g = PlasmoGraph()
g.solver = GurobiSolver(OutputFlag=0)

include("modelgen4.jl")
psize=6
otime=1:psize
oproducts=1:psize

node = []

for t in otime
    n = add_node(g)
#por que esta linea es distinta??
    m = spmodel(oproducts,t)
    setmodel(n,m)
    push!(node,n)
end
#Without Aggregation
@linkconstraint(g, [s in sites, i in oproducts, t in 1:otime[end-1]], node[t][:vf][s,i,t]  == node[t+1][:vi][s,i,t+1])
#Aggregation of multipliers - sites
#@linkconstraint(g, [i in oproducts, t in 1:otime[end-1]], sum(node[t][:vf][s,i,t] for s in sites)  == sum(node[t+1][:vi][s,i,t+1] for s in sites))
#Aggregation of multipliers - Product
#@linkconstraint(g, [s in sites, t in 1:otime[end-1]], sum(node[t][:vf][s,i,t] for i in oproducts) == sum(node[t+1][:vi][s,i,t+1] for i in oproducts))
#Aggregation of multipliers - sites & Product
#@linkconstraint(g, [t in 1:otime[end-1]], sum(node[t][:vf][s,i,t] for i in oproducts, s in sites) == sum(node[t+1][:vi][s,i,t+1] for i in oproducts, s in sites))
function cheat6(mf)
  return 72827.587
end
