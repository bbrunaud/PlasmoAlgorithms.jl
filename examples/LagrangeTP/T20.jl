using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi
using DataFrames

Logging.configure(level=DEBUG)

g = PlasmoGraph()
g.solver = GurobiSolver(OutputFlag=0)

include("modelgen4.jl")
psize=20
otime=1:psize
oproducts=1:psize
node = Dict()
for t in otime
  n = add_node(g)
  m = spmodel(oproducts,t)
  setmodel(n,m)
  node[t] = n
end


@linkconstraint(g, [s in sites, i in oproducts, t in 1:otime[end-1]],node[t][:vf][s,i,t] == node[t+1][:vi][s,i,t+1])

function cheat6(mf)
  return 72827.587
end

function heur(mf)
  return 515551.12
end
