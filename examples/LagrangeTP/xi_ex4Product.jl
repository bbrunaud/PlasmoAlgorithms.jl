G = []
#include("ex4Product.jl")
#push!(G,g)
include("ex4Temporal.jl")
push!(G,g)
#include("ex4ProdTemp.jl")
#push!(G,g)

names = ["Temporal","ProdTemp"]

df = DataFrame(Xi1=[],Xi2=[],Gap=[],LastIter=[],Time=[],Problem=[])

x1 = 0.05
x2 = 0.25

for i in 1:length(G)
  for x1 in 0.025:0.025:0.1
    for x2 in 0.5:0.5
      res, = lagrangesolve(deepcopy(G[i]),ξ1=x1,ξ2=x2,max_iterations=50,solveheuristic=cheat6)
      push!(df,[x1,x2,res[:Gap],res[:Iterations],res[:Time],names[i]])
      writetable("xi_Temporal.csv",df)
    end
  end
end


