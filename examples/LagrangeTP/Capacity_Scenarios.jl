using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi
using DataFrames
using JLD

Logging.configure(level=DEBUG)



psize=20
timeperiods=1:psize
products=1:psize



iters=[]
gaps =[]
solvetimes =[]
objvals =[]
bestbounds=[]

heuristics=[]
for k in 1:10
#k=10
    g = PlasmoGraph()
    g.solver = GurobiSolver(OutputFlag=0)
    include("modelgen4.jl")
    HT=720*(1+k/100)
    #Create Monolitic Model
    m=spmodel(products,timeperiods)
    #Solve model
    solve(m)

    Heuristic=getobjectivevalue(m)
    #558486.802595468
    #
    println(Heuristic)
    node = []
        for i in products
            n = add_node(g)
            m = spmodel(i)
            setmodel(n,m)
            push!(node,n)
        end

    @linkconstraint(g, [s in sites, i in 1:products[end-1], t in timeperiods],
    node[i][:hf][s,i,t] == node[i+1][:hi][s,i+1,t])


    if 6 in psize
        function heur(mf)
            return 72827.587
        end
    elseif 20 in psize
        function heur(mf)
            return Heuristic
        end
    end

    method = :subgradient
    δ = 0.9
    maxiter = 50
    timelimit = 1000

    result = lagrangesolve(deepcopy(g),max_iterations=maxiter,update_method=method,δ=δ,timelimit=timelimit, lagrangeheuristic=heur, initialmultipliers=:relaxation)
    #,solveheuristic=cheat20, λinit=:zero)

    #iters = result.numiterations
    push!(iters,result.numiterations)
    push!(gaps,result.gap)
    push!(solvetimes,result.solvetime)
    push!(objvals,result.objval)
    push!(bestbounds, result.bestbound)
    push!(heuristics,Heuristic)
    #save("gaps",gaps)
    #save("solvetimes", solvetime)
    #save("objvals", objval)
end
