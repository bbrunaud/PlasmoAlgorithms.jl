using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi

Logging.configure(level=DEBUG)

g = PlasmoGraph()
g.solver = GurobiSolver(OutputFlag=0)

include("modelgen.jl")
node = Dict()
for t in otime
  n = add_node(g)
  m = spmodel(oproducts,t)
  setmodel(n,m)
  node[t] = n
end


@linkconstraint(g, [s in sites, i in oproducts, t in 1:otime[end-1]],node[t][:vf][s,i,t] == node[t+1][:vi][s,i,t+1])
