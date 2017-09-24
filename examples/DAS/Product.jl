using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi

Logging.configure(level=DEBUG)

g = PlasmoGraph()
g.solver = GurobiSolver(MIPGap=0.01)

include("ProductDecomposition.jl")

customers = ["CUS$i" for i in 1:10]
numproducts = 10
products = ["PROD$i" for i in 1:numproducts]
periods = 1:12

numproducts = length(products)
origins = vcat(plants,warehouses)
destinations = vcat(warehouses,customers)

node = []
for i in 1:numproducts
    n = add_node(g)
    m = model(customers,products[i],1:12)
    setmodel(n,m)
    push!(node,n)
end

# Linkings
@linkconstraint(g, [i in plants, p in 1:(numproducts-1), t in periods], node[p][:capf][i,products[p],t] == node[p+1][:capi][i,products[p+1],t])
@linkconstraint(g, [j in warehouses, p in 1:(numproducts-1), t in periods], node[p][:y][j,t,products[p]] == node[p+1][:y][j,t,products[p+1]])
@linkconstraint(g, [i in origins, j in destinations, p in 1:(numproducts-1), t in periods], node[p][:unf][i,j,products[p],t] == node[p+1][:uni][i,j,products[p+1],t])

function fixwarehouses(mf)
  for j in warehouses
    for t in perioods
      vars = []
      vals = []
      for p in 1:numproducts
        push!(vars,getindex(mf,Symbol("node$p.y[$j,$t,PROD$p]")))
        push!(vals,getvalue(vars[p]))
      end
      maxy = maximum(vals)
      map(setlowerbound,vars,[maxy for k in 1:length(vars)])
    end
  end
  status = solve(mf)
  if status == :Optimal
    return getobjectivevalue(mf)
  else
    error("Heuristic model not infeasible or unbounded")
  end
end

result = lagrangesolve(g,solveheuristic=fixwarehouses)

display(result)
println("")
