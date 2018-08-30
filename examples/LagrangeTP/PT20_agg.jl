using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi
using DataFrames

Logging.configure(level=DEBUG)

g = PlasmoGraph()
g.solver = GurobiSolver(OutputFlag=0)

include("modelgen4.jl")
#psize=6
psize=20
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



#Without Aggregation
@linkconstraint(g, [s in sites, i in 1:oproducts[end-1], t in otime],node[i,t][:hf][s,i,t] == node[i+1,t][:hi][s,i+1,t])
@linkconstraint(g, [s in sites, i in oproducts, t in 1:otime[end-1]],node[i,t][:vf][s,i,t] == node[i,t+1][:vi][s,i,t+1])

#=
#Aggregation multipliers Sites
for i in oproducts
  for t in otime
    m = getmodel(node[i,t])
    hi = getindex(m,:hi)
    hf = getindex(m,:hf)
    vi = getindex(m,:vi)
    vf = getindex(m,:vf)
    @variable(m, svi[p in [i],τ in [t]] >= 0)
    @variable(m, svf[p in [i],τ in [t]] >= 0)
    @constraint(m, agg1, svi[i,t] == sum(vi[s,i,t] for s in sites))
    @constraint(m, agg2, svf[i,t] == sum(vf[s,i,t] for s in sites))
    @variable(m, shi[p in [i],τ in [t]] >= 0)
    @variable(m, shf[p in [i],τ in [t]] >= 0)
    @constraint(m, agg3, shi[i,t] == sum(hi[s,i,t] for s in sites))
    @constraint(m, agg4, shf[i,t] == sum(hf[s,i,t] for s in sites))
  end
end
@linkconstraint(g, [i in 1:oproducts[end-1], t in otime], node[i,t][:shf][i,t] == node[i+1,t][:shi][i+1,t])
@linkconstraint(g, [i in oproducts, t in 1:otime[end-1]], node[i,t][:svf][i,t] == node[i,t+1][:svi][i,t+1])
=#
#function heur(mf)
#  return 72669.0622
#end

function heur(mf)
      return 515551.12
 end
