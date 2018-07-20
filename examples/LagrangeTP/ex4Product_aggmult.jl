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

for i in oproducts
    n = add_node(g)
    m = spmodel(i)
    setmodel(n,m)
    push!(node,n)
end

G = []

#Aggregation of multipliers - sites
#g1 = deepcopy(g)
#@linkconstraint(g, [i in 1:oproducts[end-1], t in otime], sum(node[i][:hf][s,i,t] for s in sites) == sum(node[i+1][:hi][s,i+1,t] for s in sites))
#push!(G,g1)

#Aggregation of multipliers - periods
#g2 = deepcopy(g)
#@linkconstraint(g, [i in 1:oproducts[end-1], s in sites], sum(node[i][:hf][s,i,t] for t in otime) == sum(node[i+1][:hi][s,i+1,t] for t in otime))
#push!(G,g2)

#Aggregation of multipliers - sites+periods
#g3 = deepcopy(g)
@linkconstraint(g, [i in 1:oproducts[end-1]], sum(node[i][:hf][s,i,t] for t in otime for s in sites) == sum(node[i+1][:hi][s,i+1,t] for t in otime for s in sites))
#push!(G,g3)

function heur(mf)
  return 72827.587
end
