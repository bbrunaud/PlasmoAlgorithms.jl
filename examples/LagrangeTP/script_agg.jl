using Plasmo
using PlasmoAlgorithms
using Gurobi
using JLD
using DataFrames

files = ["T6_agg.jl",
        "P6_agg.jl"]#,
#        "PT6_agg.jl"]
k = 3

methods = [:subgradient]#,
           #:bettersubgradient,
           #:fastsubgradient,
           #:intersectionstep,
           #:marchingstep]

initial = [:relaxation]
d = Dict()

maxiter=3000
timelimit = 90

for f in files
  for method in methods
    for init in initial
      include(f)
      #mf = create_flat_graph_model(g)
      #mf.solver = GurobiSolver()
      #solve(mf)
      #d[f,k,:relax] = getobjectivevalue(mf)
      r = lagrangesolve(g,max_iterations=maxiter,timelimit=timelimit,lagrangeheuristic=heur,initialmultipliers=init,update_method=method)
      d[f,k,method,init] = r
      d[f,k,:length] = length(getlinkconstraints(g))
    end
  end
end

save("Agg_$k.jld","d",d)
