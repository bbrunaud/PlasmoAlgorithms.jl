using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi
using DataFrames

Logging.configure(level=DEBUG)
srand(1)

include("modelgen4.jl")
psize=6
otime=1:psize
oproducts=1:psize

function heur(mf)
  return 72827.587
end

orders = []
gaps =[]

for k in 1:3

pord = vcat([1],shuffle(2:6))

if order in orders
    continue
end

g = PlasmoGraph()
g.solver = GurobiSolver(OutputFlag=0)

node1 = []

for i in oproducts
    n = add_node(g)
    m = spmodel(i)
    setmodel(n,m)
    push!(node1,n)
end

node = [node1[pord[i]] for i in 1:length(pord)]
push!(orders,pord)
println("Order = $pord")
@linkconstraint(g, [s in sites, i in 1:5, t in otime],node[i][:hf][s,pord[i],t] == node[i+1][:hi][s,pord[i+1],t])

r = lagrangesolve(g, max_iterations=20, lagrangeheuristic=heur, initialmultipliers=:relaxation)
push!(gaps,r.gap)
end
