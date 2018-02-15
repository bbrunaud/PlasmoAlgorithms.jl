G = []
include("ex4Product.jl")
push!(G,g)
include("ex4Temporal.jl")
push!(G,g)
include("ex4ProdTemp.jl")
push!(G,g)

names = ["Product","Temporal","ProdTemp"]

df = DataFrame(Xi1=[],Xi2=[],Gap=[],LastIter=[],Time=[],Problem=[])

x1 = 0.05
x2 = 0.25

for i in 1:length(G)
  for x1 in 0.01:0.05:0.1
    for x2 in 0.1:0.5:1
      res, = lagrangesolve(deepcopy(G[i]),ξ1=x1,ξ2=x2,max_iterations=2)
      push!(df,[x1,x2,res[:Gap],res[:Iterations],res[:Time],names[i]])
    end
  end
end

writetable("xi_Product.csv",df)
