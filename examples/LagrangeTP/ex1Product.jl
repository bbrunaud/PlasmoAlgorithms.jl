using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi

Logging.configure(level=DEBUG)

g = PlasmoGraph()
g.solver = GurobiSolver(OutputFlag=0)

include("modelgen.jl")
node = []
for i in oproducts
    n = add_node(g)
    m = spmodel(i)
    setmodel(n,m)
    push!(node,n)
end

@linkconstraint(g, [s in sites, i in 1:oproducts[end-1], t in otime],node[i][:hf][s,i,t] == node[i+1][:hi][s,i+1,t])

result, df = lagrangesolve(g)

display(result)
println("")
