# Parallel model solve function, returns an array of objective values with dimension equal to of elements in the collection for which pmap was applied
function psolve(m::JuMP.Model)
  solve(m)
  d = Dict()
  d[:objective] = getobjectivevalue(m)
  d[:values] = m.colVal
  node = getnode(m)
  for v in values(node.index)
    d[:nodeindex] = v
  end
  println("Solved node $(d[:nodeindex]) on $(gethostname())")
  return d
end

function  lagrangesolve(graph::PlasmoGraph;update_method=:subgradient_original,max_iterations=100,ϵ=0.001,α=2,UB=5e5,LB=-1e5,δ=0.95)
  tic()
  starttime = time()
  df = DataFrame(Iter=[],Time=[],α=[],step=[],UB=[],LB=[],Hk=[],Zk=[],Gap=[])
  res = Dict()
  # 1. Check for dynamic structure. If not error
  # TODO

  # 2. Generate model for heuristic
  mflat = create_flat_graph_model(graph)
  mflat.solver = graph.solver

  # 2. Generate master problem
  ## Number of multipliers
  links = getlinkconstraints(graph)
  nmult = length(links)

  ## Master Model
  ms = Model(solver=graph.solver)
  @variable(ms, η)
  @variable(ms, λ[1:nmult])
  @objective(ms, Min, η)

  # Equality constraint the multiplier is unbounded in sign. For <= or >= need to set the lower or upper bound at 0
  # TODO ... only accepting equality constraints with rhs 0 by now.

  # 3. Generate subproblem array
  # Assuming nodes are created in order
  SP = [graph.nodes[i].attributes[:model] for i in 1:length(graph.nodes)]
  for sp in SP
    JuMP.setsolver(sp,graph.solver)
  end
  SPObjectives = [graph.nodes[i].attributes[:model].obj for i in 1:length(graph.nodes)]
  sense = SP[1].objSense
  # To update the multiplier in the suproblem, call @objective again

  # 4. Initialize
####
  θ=0
####
  Kprev=[0 for j in 1:nmult]
  λk = [0 for j in 1:nmult]
  i = 0

  # 5. Solve subproblems
  for iter in 1:max_iterations
    debug("*********************")
    debug("*** ITERATION $iter  ***")
    debug("*********************")
    Zprev = 0
    SPR = pmap(psolve,SP)
    Zk = 0
    nodedict = getnodes(graph)
    for spd in SPR
      Zk += spd[:objective]
      getmodel(nodedict[spd[:nodeindex]]).colVal = spd[:values]
    end


    Zprev = Zk
    debug("Zk = $Zk")

    # 7. Solve Lagrange heuristic (fix integers and binaries)
    mflat.colVal = vcat([getmodel(n).colVal for n in values(nodedict)]...)
    for j in 1:mflat.numCols
      if mflat.colCat[j] in (:Bin,:Int)
        mflat.colUpper[j] = mflat.colVal[j]
        mflat.colLower[j] = mflat.colVal[j]
      end
    end
    solve(mflat)
    Hk = getobjectivevalue(mflat)
    # While mflat switches from Max to Min
    Hk *= sense==:Min ? 1 : -1
    debug("Hk = $Hk")

    # 8. Update bounds and check bounds convergence
    # Minimization problem
###
    UBprev = UB
###
    LBprev = LB
    if sense == :Min
      if iter > 1
        i += LB == LBprev ? 1 : -i
        α *= i>2 ? 0.95 : 1
      end
      LB = max(Zk,LB)
      UB = min(Hk,UB)
      graph.objVal = UB
      gap = (UB - LB)/UB
      gap < ϵ &&  break
    else
      if iter > 1
        i += UB == UBprev ? 1 : -i
        α *= i>2 ? δ : 1
      end
      LB = max(Hk,LB)
      UB = min(Zk,UB)
      graph.objVal = LB
      gap = (UB - LB)/LB
      gap < ϵ &&  break
    end


    #if α < 0.001 break
    #end
    res[:Objective] = sense == :Min ? UB : LB
    res[:BestBound] = sense == :Min ? LB : UB
    res[:Iterations] = iter


    # 9. Update λ
    λprev = λk

    # 9. Multipliers Update
    lval = [getvalue(links[j].terms) for j in 1:nmult]
    dif = sense == :Max ? Zk-LB : UB - Zk
    μ=lval+θ* Kprev
#    if norm(Kprev)>0 && dot(lval,Kprev)<0
    if  dot(lval,Kprev)<0
        θ=norm(lval)/norm(Kprev)
      else
        θ=0
      end



    # update multiplier bounds (Bundle method)
    if update_method == :bundle
      for j in 1:nmult
        setupperbound(λ[j], λprev[j] + step*abs(lval[j]))
        setlowerbound(λ[j], λprev[j] - step*abs(lval[j]))
      end
    end
    # Cutting planes or Bundle
    if update_method in (:cuttingplanes,:bundle)
      @constraint(ms, η >= Zk + sum(λ[j]*lval[j] for j in 1:nmult))
      debug("Last cut = $(ms.linconstr[end])")
      solve(ms)
      λk = getvalue(λ)
    end
    # Subgradient
    if update_method == :subgradient
      step = α*dif/dot(lval,μ)
      λk = λprev - step*μ
    end

    if update_method == :subgradient_original
      step = α*dif/norm(lval)^2
      λk = λprev - step*lval
    end

    debug("Step = $step")
    debug("α = $α")
    # 10. Update objectives
    # Restore initial objective
    for (j,sp) in enumerate(SP)
      sp.obj = SPObjectives[j]
    end
    # add dualized part
    for l in 1:nmult
      for j in 1:length(links[l].terms.vars)
        var = links[l].terms.vars[j]
        coeff = links[l].terms.coeffs[j]
        var.m.obj += λk[l]*coeff*var
      end
    end

    debug("UB = $UB")
    debug("LB = $LB")
    debug("gap = $gap")
    push!(df,[iter,round(time()-starttime),α,step,UB,LB,Hk,Zk,gap])
  end
  res[:Time] = toc()
  return res, df
end
