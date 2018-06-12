include("ex4Product.jl")

method = :subgradient
δ = 0.9
maxiter = 100
timelimit = 1000



result, df = lagrangesolve(deepcopy(g),max_iterations=maxiter,update_method=method,δ=δ,timelimit=timelimit, lagrangeheuristic=heur, initialmultipliers=:relaxation)
#,solveheuristic=cheat20, λinit=:zero)

iters = result[:Iterations]
df[:Example] = ["Example4 Product" for i in 1:iters]
df[:Method] = [method for i in 1:iters]
df[:δ] = [δ for i in 1:iters]

result[:Example] = "Example4 Product"
result[:Method] = method
result[:δ] = δ

writetable("res20zero_ex4Product.csv",df)
