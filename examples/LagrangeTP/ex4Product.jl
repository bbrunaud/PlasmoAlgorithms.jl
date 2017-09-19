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

for i in oproducts
    n = add_node(g)
    m = spmodel(i)
    setmodel(n,m)
    push!(node,n)
end

@linkconstraint(g, [s in sites, i in 1:oproducts[end-1], t in otime],node[i][:hf][s,i,t] == node[i+1][:hi][s,i+1,t])

methods = [:subgradient_original,:subgradient]
#methods = [:subgradient_original]
Δ = 0.5:0.1:1.0
#Δ = 0.95
maxiter = 1000

DF = DataFrame(Iter=[],Time=[],α=[],step=[],UB=[],LB=[],Hk=[],Zk=[],Gap=[],Example=[],Method=[],δ=[])
for method in methods
  for δ in Δ
    result, df = lagrangesolve(deepcopy(g),max_iterations=maxiter,update_method=method,δ=δ)
    iters = result[:Iterations]
    df[:Example] = ["Example4 Product" for i in 1:iters]
    df[:Method] = [method for i in 1:iters]
    df[:δ] = [δ for i in 1:iters]
    DF = vcat(DF,df)
  end
end


writetable("ex4Product.csv", DF)
