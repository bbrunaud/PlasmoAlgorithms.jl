using Plasmo
using PlasmoAlgorithms
using Gurobi
using JLD
using DataFrames

n = 6

files = ["T$n.jl",
        "P$n.jl",
        "PT$n.jl"]

methods = [:subgradient,
           :bettersubgradient,
           :fastsubgradient,
           :intersectionstep,
           :marchingstep]

initial = [:zero, :relaxation]
d = Dict()

maxiter=30

for f in files
  for method in methods
    for init in initial
      include(f)
      r = lagrangesolve(g,max_iterations=maxiter,lagrangeheuristic=heur,initialmultipliers=init,update_method=method)
      d[f,method,init] = r
    end
  end
end

save("methods_$n.jld","d",d)
