using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi
using DataFrames

Logging.configure(level=DEBUG)

g = PlasmoGraph()
g.solver = GurobiSolver(OutputFlag=0)

include("modelgen4.jl")
#psize=6
psize=6
otime=1:psize
oproducts=1:psize

node = []

for i in oproducts
    n = add_node(g)
    m = spmodel(i,otime)
    setmodel(n,m)
    push!(node,n)
end
#without Aggregation
#@linkconstraint(g, [s in sites, i in 1:oproducts[end-1], t in otime], node[i][:hf][s,i,t]  == node[i+1][:hi][s,i+1,t])
#Aggregation of multipliers - sites
#@linkconstraint(g, [i in 1:oproducts[end-1], t in otime], sum(node[i][:hf][s,i,t] for s in sites) == sum(node[i+1][:hi][s,i+1,t] for s in sites))

#Aggregation of multipliers - time
#@linkconstraint(g, [s in sites, i in 1:oproducts[end-1]], sum(node[i][:hf][s,i,t] for t in otime) == sum(node[i+1][:hi][s,i+1,t] for t in otime))

#Aggregation of multipliers - sites & time
@linkconstraint(g, [i in 1:oproducts[end-1]], sum(node[i][:hf][s,i,t] for s in sites, t in otime) == sum(node[i+1][:hi][s,i+1,t] for s in sites, t in otime))

#function cheat6(mf)
#  return 72827.587
#end
