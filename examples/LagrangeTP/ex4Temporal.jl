using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi

Logging.configure(level=DEBUG)

g = PlasmoGraph()
g.solver = GurobiSolver(OutputFlag=0)

include("modelgen4.jl")
node = Dict()
for t in otime
  n = add_node(g)
  m = spmodel(oproducts,t)
  setmodel(n,m)
  node[t] = n
end


@linkconstraint(g, [s in sites, i in oproducts, t in 1:otime[end-1]],node[t][:vf][s,i,t] == node[t+1][:vi][s,i,t+1])

methods = [:subgradient_original,:subgradient]
Δ = 0.5:0.1:1.0
maxiter = 500

DF = DataFrame(Iter=[],Time=[],α=[],step=[],UB=[],LB=[],Hk=[],Zk=[],Gap=[],Example=[],Method=[],δ=[])
for method in methods
  for δ in Δ
    result, df = lagrangesolve(g,max_iterations=maxiter,update_method=method,δ=δ)
    df[:Example] = ["Example4 Temporal" for i in 1:maxiter]
    df[:Method] = [method for i in 1:maxiter]
    df[:δ] = [δ for i in 1:maxiter]
    DF = vcat(DF,df)
  end
end


writetable("ex4Temporal.csv", DF)
