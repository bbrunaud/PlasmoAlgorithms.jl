include("ex4Temporal.jl")

method = :subgradient
δ = 0.9
maxiter = 30

result, df = lagrangesolve(deepcopy(g),max_iterations=maxiter,update_method=method,δ=δ)
iters = result[:Iterations]
df[:Example] = ["Example4 Temporal" for i in 1:iters]
df[:Method] = [method for i in 1:iters]
df[:δ] = [δ for i in 1:iters]

result[:Example] = "Example4 Temporal"
result[:Method] = method
result[:δ] = δ

result
