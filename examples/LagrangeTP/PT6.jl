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

@linkconstraint(g, [s in sites, i in 1:oproducts[end-1], t in otime],node[i,t][:hf][s,i,t] == node[i+1,t][:hi][s,i+1,t])
@linkconstraint(g, [s in sites, i in oproducts, t in 1:otime[end-1]],node[i,t][:vf][s,i,t] == node[i,t+1][:vi][s,i,t+1])

function heur(mf)
  return 72827.587
end

function cheat20(mf)
  return 515551.12
end
