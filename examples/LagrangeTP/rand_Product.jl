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

df = DataFrame(Order=[],Gap=[],LastIter=[],Time=[])

maxiter = 2
for k in 1:maxiter
  A = collect(oproducts)
  B = A[randperm(length(A))]

  gi = deepcopy(g)
  nodes = getnodes(gi)
  @linkconstraint(gi, [s in sites, i in 1:(length(B)-1), t in otime],nodes[B[i]][:hf][s,B[i],t] == nodes[B[i+1]][:hi][s,B[i+1],t])
  res, = lagrangesolve(g,max_iterations=4)
  push!(df,[B,res[:Gap],res[:Iterations],res[:Time]])
end

function cheat6(mf)
  return 72827.587
end
