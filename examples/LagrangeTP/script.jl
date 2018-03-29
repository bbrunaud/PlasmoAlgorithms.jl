using Plasmo
using PlasmoAlgorithms
using Gurobi
using JLD
using DataFrames

n = 3

files = ["T$n.jl",
        "P$n.jl",
        "PT$n.jl"]

methods = [:subgradient]#,
           #:bettersubgradient,
           #:fastsubgradient,
           #:intersectionstep,
           #:marchingstep]

initial = [:zero, :relaxation]
d = Dict()

maxiter=3000
timelimit = 2

for f in files
  for method in methods
    for init in initial
      include(f)
      r = lagrangesolve(g,max_iterations=maxiter,timelimit=timelimit,lagrangeheuristic=heur,initialmultipliers=init,update_method=method)
      d[f,method,init] = r
    end
  end
end

save("Standard_$n.jld","d",d)
