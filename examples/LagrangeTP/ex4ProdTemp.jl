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

methods = [:subgradient_original,:subgradient]
Δ = 0.5:0.1:1.0
maxiter = 1000

DF = DataFrame(Iter=[],Time=[],α=[],step=[],UB=[],LB=[],Hk=[],Zk=[],Gap=[],Example=[],Method=[],δ=[])
for method in methods
    for δ in Δ
  result, df = lagrangesolve(deepcopy(g),max_iterations=maxiter,update_method=method,δ=δ)
  iters = result[:Iterations]
  df[:Example] = ["Example4 Product-Temporal" for i in 1:iters]
  df[:Method] = [method for i in 1:iters]
  df[:δ] = [δ for i in 1:iters]
  DF = vcat(DF,df)
  end
end

writetable("ex4ProdTemp.csv", DF)
