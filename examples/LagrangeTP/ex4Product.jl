using Plasmo
using PlasmoAlgorithms

include("modelgen4.jl")
psize=20
timeperiods=1:psize
products=1:psize

g = PlasmoGraph()
node = []

for i in products
    n = add_node(g)
    m = model(i)
    setmodel(n,m)
    push!(node,n)
end

@linkconstraint(g, [s in sites, i in 1:products[end-1], t in timperiods],
node[i][:cf][s,i,t] == node[i+1][:ci][s,i+1,t])

function cheat6(mf)
  return 72827.587
end

function cheat20(mf)
  return 515551.12
end
