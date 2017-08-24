using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi

Logging.configure(level=DEBUG)

g = PlasmoGraph()
g.solver = GurobiSolver(OutputFlag=0)

include("modelgen4.jl")
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

result = lagrangesolve(g)

display(result)
println("")
