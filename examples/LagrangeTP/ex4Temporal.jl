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
for t in otime
  n = add_node(g)
  m = spmodel(oproducts,t)
  setmodel(n,m)
  node[t] = n
end


@linkconstraint(g, [s in sites, i in oproducts, t in 1:otime[end-1]],node[t][:vf][s,i,t] == node[t+1][:vi][s,i,t+1])

methods = [:subgradient]
Δ = [0.5]
maxiter = 50

DF = DataFrame(Iter=[],Time=[],α=[],step=[],UB=[],LB=[],Hk=[],Zk=[],Gap=[],Example=[],Method=[],δ=[])
for method in methods
  for δ in Δ
    result, df = lagrangesolve(g,max_iterations=maxiter,update_method=method,δ=δ)
    iters = result[:Iterations]
    df[:Example] = ["Example4 Temporal" for i in 1:iters]
    df[:Method] = [method for i in 1:iters]
    df[:δ] = [δ for i in 1:iters]
    DF = vcat(DF,df)
  end
end


writetable("ex4Temporal.csv", DF)
