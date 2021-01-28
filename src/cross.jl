
function crossprepare(c::CrossGraph, cuts=[:LP])
    # cross-specific mappings
    # Lagrange mappings in Benders Graph
    links = getlinkconstraints(c.bd)
    for n in getnodes(c.bd)
        mn = getmodel(n)
        mn.ext[:preobj] = objective_function(mn)
        mn.ext[:multmap] = Dict()
        mn.ext[:varmap] = Dict()
    end
    for (i,link) in enumerate(links)
        for (j,term) in enumerate(link.func.terms)
            var = term[1]
            coeff = term[2]
            owner_model(var).ext[:multmap][i] = (coeff, var)
            owner_model(var).ext[:varmap][var] = (i,j)
          end
    end

    # Benders mappings in Lagrange Graph
    PlasmoAlgorithms.identifylevels(c.lg)
    links = getlinkconstraints(c.lg)
    for (numlink,link) in enumerate(links)
        #Take the two variables of the constraint
        var1, var2 = collect(keys(link.func.terms))
        #Determine which nodes they belong to
        nodeV1 = getnode(var1)
        nodeV2 = getnode(var2)
        #Set the order of the nodes
        if nodeV1 in out_neighbors(c.lg,nodeV2)
            childnode = nodeV1
            childvar = var1
            parentnode = nodeV2
            parentvar = var2
        else
            childnode = nodeV2
            childvar = var2
            parentnode = nodeV1
            parentvar = var1
        end
        childindex = c.lg[childnode]
        push!(parentnode[:childvars][childindex], parentvar)
        if !hasattribute(parentnode, :λcomponents)
            parentnode[:λcomponents] =  Dict()
        end
        if !haskey(parentnode[:λcomponents], childindex)
            parentnode[:λcomponents][childindex] = Int64[]
        end
        push!(parentnode[:λcomponents][childindex], numlink)
    end

    # Standard preparations
    PlasmoAlgorithms.bdprepare(c.bd)
    PlasmoAlgorithms.lgprepare(c.lg, δ=0.5, maxnoimprove=3,cpbound=1e6)
 end


 function crossoptimize!(g::OptiGraph; max_iterations=10, subgradientiterations=5, ϵ=0.005, timelimit=200, bdcuts=[:LP], singleline=true)
     c = CrossGraph(g, deepcopy(g))
     crossoptimize!(c, max_iterations=max_iterations, subgradientiterations=subgradientiterations, ϵ=ϵ, timelimit=timelimit, bdcuts=bdcuts, singleline=singleline)
 end


 function crossoptimize!(c::CrossGraph; max_iterations=10, subgradientiterations=5, ϵ=0.005, timelimit=200, bdcuts=[:LP],singleline=true, verbose=true)
     s = Solution(method=:cross)
     starttime = time()
     crossprepare(c, bdcuts)
     n = getattribute(c.bd, :normalized)

     for i in 1:max_iterations
         tic = time()
         bs = bendersoptimize!(c.bd, max_iterations=1, cuts=bdcuts, verbose=false)
         benders_to_lagrange(c)
         updatemethod =  i <= subgradientiterations ? :subgradient : :cuttingplanes
         ls = lagrangeoptimize!(c.lg,max_iterations=1, update_method=updatemethod, lagrangeheuristic=x->bs.objval, verbose=false)
         lagrange_to_benders(c, termination=ls.termination)

         UB = bs.objval
         LB = ls.bestbound

         tstamp = time() - starttime

         itertime = time() - tic

         if n == 1
             saveiteration(s,tstamp,[UB,LB,itertime,tstamp],n)
         else
             saveiteration(s,tstamp,[n*LB,n*UB,itertime,tstamp],n)
         end
         if verbose 
            i == 1 && singleline && printheader(c.bd, :Cross)
            printiterationsummary(s,singleline=singleline)
          end

         # Check Optimality
         if abs(UB-LB)/(1e-10 + abs(UB)) < ϵ
             s.termination = "Optimal"
             return s
         end

         # Check time limit
         if tstamp > timelimit
           s.termination = "Time Limit"
           return s
         end
     end
     s.termination = "Max Iterations"
     return s
 end

function benders_to_lagrange(c::CrossGraph)
    # Push Cut Data
    cutdataarray = getattribute(c.lg, :cutdata)
    for node in getnodes(c.bd)
        nodeindex = c.bd[node]
        m = getmodel(node)
        preobjval = value(m.ext[:preobj])
        λcomponent = Int64[]
        coeffs = Float64[]
        xk = Float64[]
        for k in keys(m.ext[:multmap])
            coeff = m.ext[:multmap][k][1]
            var = m.ext[:multmap][k][2]
            push!(λcomponent, k)
            push!(coeffs, coeff)
            push!(xk, JuMP.getvalue(var))
        end
        push!(cutdataarray[nodeindex], LagrangeCutData(preobjval, λcomponent, coeffs, xk))
    end
end

function lagrange_to_benders(c::CrossGraph; termination="Max Iterations")
    if termination == "Max Iterations"
        λ = getattribute(c.lg,:λ)[end-1]
        println("Max Iterations $λ")
    else
        λ = getattribute(c.lg,:λ)[end]
        println("Optimal $λ")
    end
    for lgnode in getnodes(c.lg)
        if in_degree(c.lg, lgnode) > 0
            parentnode = in_neighbors(c.lg, lgnode)[1]
            nodeindex = c.lg[lgnode]
            parentindex = c.lg[parentnode]
            λk = λ[parentnode[:λcomponents][nodeindex]]
            bdnode = getnode(c.bd, parentindex)
            zld = PlasmoAlgorithms.subgraphobjective(lgnode, c.lg)
            # TODO: Correct Cut Data type and push
            # ERROR: MethodError: no method matching push!(::Distributed.RemoteChannel{Channel{PlasmoAlgorithms.CutData}}, ::PlasmoAlgorithms.LagrangeCrossCutData)
            put!(bdnode[:cutdata][nodeindex], LagrangeCrossCutData(zld,λk))
        end
    end
end

