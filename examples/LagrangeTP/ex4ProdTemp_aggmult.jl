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

node = Dict()
for i in oproducts
  for t in otime
    n = add_node(g)
    m = spmodel(i,t)
    setmodel(n,m)
    node[i,t] = n
  end
end

@linkconstraint(g, [i in 1:oproducts[end-1], t in otime], sum(node[i,t][:hf][s,i,t] for s in sites) == sum(node[i+1,t][:hi][s,i+1,t] for s in sites))
@linkconstraint(g, [s in sites, t in 1:otime[end-1]], sum(node[i,t][:vf][s,i,t] for i in oproducts) == sum(node[i,t+1][:vi][s,i,t+1] for i in oproducts))
