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
K=2
for i in 1:K:length(oproducts)
    n = add_node(g)
    m = spmodel(i:(i+1),otime)
    setmodel(n,m)
    push!(node,n)
end

@linkconstraint(g, [s in sites, i in 1:(length(node)-1), t in otime],node[i][:hf][s,K*i,t] == node[i+1][:hi][s,K*i+1,t])

function cheat6(mf)
  return 72827.587
end
