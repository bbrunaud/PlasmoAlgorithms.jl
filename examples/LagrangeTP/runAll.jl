using DataFrames

filenames = ["ex4Product", "ex4Temporal", "ex4ProdTemp"]

method = :subgradient
δf = 0.9
maxiter = 300
timelimit = 200




out = DataFrame(Example=[],Size=[],init=[],LB=[],UB=[],Gap=[],Time=[],Iterations=[])
for f in filenames
    for s in (6,20)
        for linit = (:zero,:relaxation)

            include("$f.jl")

            cheatfun(mf)  = eval(Expr(:call,Symbol("cheat$s"),mf))


            g = graphgen(s)

            result, df = lagrangesolve(deepcopy(g),max_iterations=maxiter,update_method=method,δ=δf,timelimit=timelimit,solveheuristic=cheatfun, λinit=linit)
            iters = result[:Iterations]
            df[:Example] = [f for i in 1:iters]
            df[:Method] = [method for i in 1:iters]
            df[:δ] = [δ for i in 1:iters]

            result[:Example] = f
            result[:Method] = method
            result[:δ] = δ

            push!(out,[f,s,linit,result[:BestBound],result[:Objective],result[:Gap],result[:Time],result[:Iterations]])
            show(out)

        end
    end
end

writetable("res.csv",out)
