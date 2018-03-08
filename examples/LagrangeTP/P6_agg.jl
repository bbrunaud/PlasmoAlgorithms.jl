using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi
using JuMP
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


#=
#Aggregation of multipliers - sites
for i in oproducts
    m = getmodel(node[i])
    hi = getindex(m,:hi)
    hf = getindex(m,:hf)
    @variable(m, shi[i in oproducts, t in otime] >= 0)
    @variable(m, shf[i in oproducts, t in otime] >= 0)
    @constraint(m, agg1[t in otime], shi[i,t] == sum(hi[s,i,t] for s in sites))
    @constraint(m, agg2[t in otime], shf[i,t] == sum(hf[s,i,t] for s in sites))
end
@linkconstraint(g, [i in 1:oproducts[end-1], t in otime], node[i][:shf][i,t] == node[i+1][:shi][i+1,t])


#Aggregation of multipliers - periods
for i in oproducts
    m = getmodel(node[i])
    hi = getindex(m,:hi)
    hf = getindex(m,:hf)
    @variable(m, shi[s in sites, i in oproducts] >= 0)
    @variable(m, shf[s in sites, i in oproducts] >= 0)
    @constraint(m, agg1[s in sites], shi[s,i] == sum(hi[s,i,t] for t in otime))
    @constraint(m, agg2[s in sites], shf[s,i] == sum(hf[s,i,t] for t in otime))
end
@linkconstraint(g, [i in 1:oproducts[end-1], s in sites], node[i][:shf][s,i] == node[i+1][:shi][s,i+1])
=#

#Aggregation of multipliers - sites+periods
for i in oproducts
    m = getmodel(node[i])
    hi = getindex(m,:hi)
    hf = getindex(m,:hf)
    @variable(m, shi[i in oproducts] >= 0)
    @variable(m, shf[i in oproducts] >= 0)
    @constraint(m, agg1, shi[i] == sum(hi[s,i,t] for t in otime for s in sites))
    @constraint(m, agg2, shf[i] == sum(hf[s,i,t] for t in otime for s in sites))
end
@linkconstraint(g, [i in 1:oproducts[end-1]], node[i][:shf][i] == node[i+1][:shi][i+1])
#=
=#

function heur(mf)
  return 72669.0622
end
