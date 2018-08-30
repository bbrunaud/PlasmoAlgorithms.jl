function solverootrelaxation(node::PlasmoNode)
  sp = getmodel(node)
  if length(sp.linconstrDuals) == 0
    solve(sp, relaxation=true)
  end
  lpfile = joinpath(tmpdir,"nodemodel.lp")
  writeLP(sp,lpfile)
  run(`cpxgetroot $lpfile 0 1`)
  lp = Model(solver=CPLEX.CplexSolver(CPX_PARAM_PREIND=0))
  lp.internalModel = MathProgBase.LinearQuadraticModel(lp.solver)
  if isfile("node0.lp")
      MathProgBase.loadproblem!(lp.internalModel,"node0.lp")
  else
    warn("Node file not found, falling back to lp relaxation")
    return solvelprelaxation(node)
  end

  # Restore Bounds
  MathProgBase.setvarLB!(lp.internalModel,sp.colLower)
  MathProgBase.setvarUB!(lp.internalModel,sp.colUpper)

  MathProgBase.optimize!(lp.internalModel)

  run(`mv node0.lp $tmpdir/`)

  dualconstraints = node.attributes[:linkconstraints]

  rootduals = MathProgBase.getconstrduals(lp.internalModel)
  sp.linconstrDuals = MathProgBase.getconstrduals(lp.internalModel)[1:length(sp.linconstrDuals)]

  λnode = getdual(dualconstraints)
  nodebound = MathProgBase.getobjval(lp.internalModel)

  node.attributes[:bound] = nodebound
  node.attributes[:λ] = λnode
end

function solverootrelaxation_debug(node::PlasmoNode, graph)
  sp = getmodel(node)
  if length(sp.linconstrDuals) == 0
    solve(sp, relaxation=true)
  end
  lpfile = joinpath(tmpdir,"nodemodel.lp")
  writeLP(sp,lpfile)
  run(`cpxgetroot $lpfile 0 1`)
  lp = Model(solver=CPLEX.CplexSolver(CPX_PARAM_PREIND=0))
  lp.internalModel = MathProgBase.LinearQuadraticModel(lp.solver)
  if isfile("node0.lp") && node.label != :ROOT
  	  # temp_s = string(0)
  	  # temp_command = string("cp node0.lp nodehistory/", temp_s)
  	  # run(`$temp_command`)
      temp = open("/tmp/RootNode/nodemodel.lp")
      original_obj = readlines(temp)[1:2]
      close(temp)
      temp = open("node0.lp")
      node0lp = readlines(temp)[3:end]
      close(temp)
      temp = open("node0.lp", "w")
      for line in original_obj
        write(temp, line)
        write(temp, "\n")
      end
      start = false
      for line in node0lp
        if line == "Subject To"
          start = true 
        end
        if start
          write(temp, line)
          write(temp, "\n")
        end
      end
      close(temp)

      MathProgBase.loadproblem!(lp.internalModel,"node0.lp")

      #check if the cuts are valid 
      all_constraints = MathProgBase.getconstrmatrix(lp.internalModel)
      all_LB =  MathProgBase.getconstrLB(lp.internalModel)
      all_UB =  MathProgBase.getconstrUB(lp.internalModel)
      optimal_colVal = graph.attributes[:optimal_solutions][node.label]
      if all_constraints.n > 135 && all_constraints.m >6
        coefficients = zeros(135)
        for j in 1:135
          coefficients[j] = all_constraints[6, j]
        end
        println("coefficients of c6")
        println(coefficients)
      end
      cut_length = all_constraints.m - 35
      if cut_length >0
        for i in 1:cut_length
          expr = 0
          coefficients = zeros(135)
          for j in 1:135
            expr += optimal_colVal[j] * all_constraints[i+35, j]
            coefficients[j] = all_constraints[i+35, j]
          end
          if expr < all_LB[i+35] || expr > all_UB[i+35]
            cp("node0.lp", "nodehistory/node0.lp", remove_destination=true)
            cp("/tmp/RootNode/nodemodel.lp", "nodehistory/nodemodel.lp", remove_destination=true)
            println("optimal value")
            println(optimal_colVal)
            println("coefficients")
            println(coefficients)
            println("rhs")
            println(expr)
            println("lhs")
            println(all_UB[i+35])
            println(node.label)
            println(i)
            # all_UB[i+35] = 1e5

            error("invalid cut")
          end
        end
      end
      # MathProgBase.setconstrUB!(lp.internalModel, all_UB)
  else
    warn("Node file not found, falling back to lp relaxation")
    return solvelprelaxation(node)
  end

  # Restore Bounds
  MathProgBase.setvarLB!(lp.internalModel,sp.colLower)
  MathProgBase.setvarUB!(lp.internalModel,sp.colUpper)

  MathProgBase.optimize!(lp.internalModel)

  run(`mv node0.lp $tmpdir/`)

  dualconstraints = node.attributes[:linkconstraints]

  rootduals = MathProgBase.getconstrduals(lp.internalModel)
  sp.linconstrDuals = MathProgBase.getconstrduals(lp.internalModel)[1:length(sp.linconstrDuals)]

  λnode = getdual(dualconstraints)
  nodebound = MathProgBase.getobjval(lp.internalModel)

  node.attributes[:bound] = nodebound
  node.attributes[:λ] = λnode
  node.attributes[:rootprob] = lp 
end

