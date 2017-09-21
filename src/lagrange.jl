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
  #println("Solved node $(d[:nodeindex]) on $(gethostname())")
  return d
end

function fixbinaries(mflat,cat=[:Bin])
  for j in 1:mflat.numCols
    if mflat.colCat[j] in cat
      mflat.colUpper[j] = mflat.colVal[j]
      mflat.colLower[j] = mflat.colVal[j]
    end
  end
  status = solve(mflat)
  if status == :Optimal
    return getobjectivevalue(mflat)
  else
    error("Heuristic model not infeasible or unbounded")
  end
end

function fixbinariesandintegers(mflat)
  fixbinaries(mflat,[:Bin,:Int])
end

function  lagrangesolve(graph::PlasmoGraph;
  update_method=:subgradient,
  max_iterations=100,
  ϵ=0.001,
  α=2,
  UB=5e5,
  LB=-1e5,
  δ=0.9,
  λinit=:relaxation,
  solveheuristic=fixbinaries)

  ########## 1. Initialize ########
  tic()
  starttime = time()
  # Results outputs
  df = DataFrame(Iter=[],Time=[],α=[],step=[],UB=[],LB=[],Hk=[],Zk=[],Gap=[])
  res = Dict()

    # Get Linkings
  links = getlinkconstraints(graph)
  nmult = length(links)


  # Equality constraint the multiplier is unbounded in sign. For <= or >= need to set the lower or upper bound at 0
  # TODO ... only accepting equality constraints with rhs 0 by now.

  # Generate subproblem array
  # Assuming nodes are created in order
  SP = [graph.nodes[i].attributes[:model] for i in 1:length(graph.nodes)]
  for sp in SP
    JuMP.setsolver(sp,graph.solver)
  end
  # Capture objectives
  SPObjectives = [graph.nodes[i].attributes[:model].obj for i in 1:length(graph.nodes)]
  sense = SP[1].objSense

  # Variables
  θ = 0
  Kprev = [0 for j in 1:nmult]
  i = 0

  # Generate model for heuristic
  mflat = create_flat_graph_model(graph)
  mflat.solver = graph.solver
  # Restore mflat sense
  if sense == :Max
    mflat.objSense = :Max
    mflat.obj = -mflat.obj
  end
  # Solve realaxation
  solve(mflat,relaxation=true)
  bestbound = getobjectivevalue(mflat)
  if sense == :Max
    UB = bestbound
  else
    LB = bestbound
  end
  debug("Solved LP relaxation with value $bestbound")
  # Set starting λ to the duals of the LP relaxation
  # TODO handle NLP relaxation
  λk = λinit == :relaxation ? mflat.linconstrDuals[end-nmult+1:end] : λk = [1.0 for j in 1:nmult]
  λprev = λk

  # Master Model (generate only for cutting planes or bundle methods)
  if update_method in [:cuttingplanes,:bundle]
    ms = Model(solver=graph.solver)
    @variable(ms, η)
    @variable(ms, λ[1:nmult])
    mssense = :Min
    if sense == :Min
      mssense = :Max
    end
    @objective(ms, mssense, η)
  end

  # 5. Solve subproblems
  for iter in 1:max_iterations
    debug("*********************")
    debug("*** ITERATION $iter  ***")
    debug("*********************")

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

    # 7. Solve Lagrange heuristic
    mflat.colVal = vcat([getmodel(n).colVal for n in values(nodedict)]...)
    Hk = solveheuristic(mflat)
    debug("Hk = $Hk")

    # 8. Update bounds and check bounds convergence
    # Minimization problem
###
    UBprev = UB
###
    LBprev = LB
    improved = true
    if sense == :Min
      if iter > 4
        if LB == LBprev
          i += 1
          improved = false
        else
          i = 0
          α = 2
        end
        if i > 2
          α *= δ
          i = 0
        end
      end
      LB = max(Zk,LB)
      UB = min(Hk,UB)
      graph.objVal = UB
      gap = (UB - LB)/UB
      gap < ϵ &&  break
    else
      if iter > 1
        if UB == UBprev
          i += 1
          improved = false
        else
          i = 0
          α = 2
        end
        if i > 4
          α *= δ
          i = 0
        end
      end
      LB = max(Hk,LB)
      UB = min(Zk,UB)
      graph.objVal = LB
      gap = (UB - LB)/LB
      gap < ϵ &&  break
    end




    #end
    res[:Objective] = sense == :Min ? UB : LB
    res[:BestBound] = sense == :Min ? LB : UB
    res[:Iterations] = iter


    # 9. Update λ
    if iter == 1
      λprev = λk
      println("Set λprev at iter 1")
    end
    if improved
      λprev = λk
      println("UPDATED λprev")
    end

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

    # Subgradient
    if update_method == :subgradient || update_method in [:cuttingplanes,:bundle]
      step = α*dif/dot(lval,μ)
      λk = λprev + step*μ
    end

    if update_method == :subgradient_original
      step = α*dif/norm(lval)^2
      λk = λprev + step*lval
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
      if sense == :Max
        @constraint(ms, η >= Zk + sum(λ[j]*lval[j] for j in 1:nmult))
      else
        @constraint(ms, η <= Zk + sum(λ[j]*lval[j] for j in 1:nmult))
      end
      debug("Last cut = $(ms.linconstr[end])")
      if iter > 10
        solve(ms)
        λk = getvalue(λ)
      end
    end


    debug("Step = $step")
    debug("α = $α")
    debug("UB = $UB")
    debug("LB = $LB")
    debug("gap = $gap")
    push!(df,[iter,round(time()-starttime),α,step,UB,LB,Hk,Zk,gap])

    α < 1e-8 && break
    step < 1e-8 && break
  end
  res[:Time] = toc()
  return res, df
end
