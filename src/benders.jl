
abstract type CutData end

struct BendersCutData <: CutData
  θk
  λk
  xk
end

struct LLIntegerCutData <: CutData
  θlb
  yk
end

struct IntegerCutData <: CutData
  yk
end

(==)(cd1::BendersCutData,cd2::BendersCutData) = (cd1.θk == cd2.θk) && (cd1.λk == cd2.λk) && (cd1.xk == cd2.xk)
(==)(cd1::LLIntegerCutData,cd2::LLIntegerCutData) = (cd1.θlb == cd2.θlb) &&  (cd1.yk == cd2.yk)
(==)(cd1::IntegerCutData,cd2::IntegerCutData) = (cd1.yk == cd2.yk)

"""
bendersolve
"""
function bendersolve(graph::Plasmo.PlasmoGraph; max_iterations::Int64=10, cuts::Array{Symbol,1}=[:LP], ϵ=1e-5,UBupdatefrequency=1,timelimit=3600,verbose=false, LB=NaN, is_nonconvex=false)
  starttime = time()
  graph.attributes[:is_nonconvex] = is_nonconvex
  global tmpdir = "/tmp/RootNode" # mktempdir()
  
  updatebound = true

  verbose && info("Preparing graph")
  bdprepare(graph, cuts)
  n = graph.attributes[:normalized]
  if isnan(LB)
    verbose && info("Solve relaxation and set LB")
    mf = graph.attributes[:mflat]
    solve(mf,relaxation=true)
    LB = getobjectivevalue(graph.attributes[:mflat])
  end
  UB = Inf

  # Set bound to root node
  rootnode = graph.attributes[:roots][1]
  rootmodel = getmodel(rootnode)
  # @constraint(rootmodel, rootmodel.obj.aff >= LB)

  # Begin iterations
  verbose && info("Begin iterations")
  for i in 1:max_iterations
    tic()
    updatebound = ((i-1) % UBupdatefrequency) == 0
    LB,UB = forwardstep(graph, cuts, updatebound)

    tstamp = time() - starttime

    itertime = toc()
    if n == 1
      saveiteration(graph.attributes[:solution],tstamp,[graph.attributes[:iterUB],LB,itertime,tstamp],n)
    else
      saveiteration(graph.attributes[:solution],tstamp,[n*LB,n*UB,itertime,tstamp],n)
    end
    printiterationsummary(graph.attributes[:solution],singleline=false)

    if abs(UB-LB) < ϵ
      graph.attributes[:solution].termination = "Optimal"
      return graph.attributes[:solution]
    end

    # Check time limit
    if tstamp > timelimit
      graph.attributes[:solution].termination = "Time Limit"
      return graph.attributes[:solution]
    end

    if graph.attributes[:stalled]
      graph.attributes[:solution].termination = "Stalled"
      return graph.attributes[:solution]
    end

    if graph.attributes[:fathomed_by_bound]
      graph.attributes[:solution].termination = "node fathomed by bound"
      return graph.attributes[:solution]
    end

    if graph.attributes[:is_infeasible]
      graph.attributes[:solution].termination = "Benders master problem is infeasible"
      return graph.attributes[:solution]
    end

  end

  graph.attributes[:solution].termination = "Max Iterations"
  return graph.attributes[:solution]
end

function forwardstep(graph::PlasmoGraph, cuts::Array{Symbol,1}, updatebound::Bool)
  levels = graph.attributes[:levels]
  numlevels = length(levels)
  for level in 1:numlevels
    nodeslevel = levels[level]
    for node in nodeslevel
      solveprimalnode(node,graph,cuts,updatebound)
    end
  end
  LB = graph.attributes[:LB]
  if updatebound
    iterUB = sum(node.attributes[:preobjval] for node in values(graph.nodes))
    graph.attributes[:iterUB] = iterUB
    UB = min(graph.attributes[:UB],iterUB)
    #update best feasible x
    if UB < graph.attributes[:UB]
      root = graph.attributes[:roots][1]
      if !haskey(graph.attributes, :best_feasible_x)
        graph.attributes[:best_feasible_x] = Dict()
      end
      for varname in keys(root.attributes[:varname_to_var])
        graph.attributes[:best_feasible_x][varname] = getvalue(root.attributes[:varname_to_var][varname])
      end
    end
    graph.attributes[:UB] = UB
  else
    UB = graph.attributes[:UB]
  end
  return LB,UB
end

function solveprimalnode(node::PlasmoNode, graph::PlasmoGraph, cuts::Array{Symbol,1}, updatebound::Bool)
  # 1. Add cuts
  generatecuts(node,graph)
  # 2. Take x
  takex(node, graph)
  # 3. solve
  if :LP in cuts && in_degree(graph,node) != 0
    solvelprelaxation(node)
  end

  if :GMI in cuts && in_degree(graph,node) != 0
    solvegmirelaxation(node)
  end
  
  if updatebound
    solvenodemodel(node,graph)
  end
  # 4. put x
  putx(node,graph)
  # 5. put cuts and nodebound
  putcutdata(node,graph,cuts)
end

function solvelprelaxation(node::PlasmoNode)
  model = getmodel(node)
  status = solve(model, relaxation = true)

  # if status != :Optimal
  #   println(node.attributes[:xin])
  #   error("lp relaxation not solved to optimality")
  # end

  dualconstraints = node.attributes[:linkconstraints]

  λnode = getdual(dualconstraints)
  nodebound = getobjectivevalue(model)

  node.attributes[:bound] = nodebound
  node.attributes[:λ] = λnode

  return status
end

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

function solvegmirelaxation(node::PlasmoNode)
  # model = node.attributes[:GMI]
  model = getmodel(node)
  status=solve(model, relaxation=true)
  temp_obj = getobjectivevalue(model)
  while true 
    cpx = model.internalModel.inner 

    cbasis = CPLEX.get_basis(cpx)[1]
    rbasis = CPLEX.get_basis(cpx)[2]
    senses = CPLEX.get_constr_senses(cpx)
    constr_matrix = CPLEX.get_constr_matrix(cpx)
    rhs = CPLEX.get_rhs(cpx)
    #identify xbar x, linking constraints
    link_constr_rows = []
    xbar_indices = []
    xbar_value = []
    xbar_ub = []
    xbar_lb = []
    x_indices = []

    for i in 1:constr_matrix.m 
      row = constr_matrix[i,:]
      is_one = false 
      is_negative_one = false 
      if length(nonzeros(row)) == 2
        for j in 1:2
          if nonzeros(row)[j] == 1
            is_one = true 
          end
          if nonzeros(row)[j] == -1 
            is_negative_one = true
          end
        end
      end
      if is_one && is_negative_one 
        push!(link_constr_rows, i)
        for j in 1:length(row)
          if row[j] == -1
            push!(x_indices, j)
            push!(xbar_lb, model.colLower[j])
            push!(xbar_ub, model.colUpper[j])
          end
          
          if row[j] == 1 
            push!(xbar_indices, j)
            push!(xbar_value, model.colVal[j])
          end
        end
      end
    end

    map_col_y = zeros(Int64, length(model.colVal)) # give a column return it's the i th y 
    j = 1
    for i in 1:length(model.colVal)
      if i in xbar_indices || i in x_indices
        continue
      end
      map_col_y[i] = j 
      j += 1
    end

    #only generate cuts when y are frac 
    is_y_frac = false 
    for col in 1:length(model.colVal)
      if col in xbar_indices || col in x_indices
        continue
      end
      if (model.colCat[col] == :Int || model.colCat[col] == :Bin ) && (model.colVal[col]%1) < 1-1e-5 && (model.colVal[col]%1) > 1e-5 
        is_y_frac = true 
      end
    end
    if !is_y_frac
      break 
    end

    #only generate GMI when xbar are all at their bounds 
    is_xbar_frac = false
    for i in 1:length(xbar_value)
      if xbar_value[i] != xbar_ub[i] && xbar_value[i] != xbar_lb[i]
        is_xbar_frac = true 
      end
    end
    if is_xbar_frac
      break 
    end

    num_struct_constr = constr_matrix.m - length(link_constr_rows)
    num_y = length(model.colVal) - 2*length(link_constr_rows)
    num_x = length(link_constr_rows)
    #get B and B_N
    B = zeros(Float64, num_struct_constr+num_y, num_struct_constr+num_y)
    B_N = zeros(Float64, num_struct_constr+num_y, constr_matrix.n+num_struct_constr+num_y)
    new_rhs = zeros(Float64, num_struct_constr+num_y)
    #identify the columns of the basis variables 
    col_cbasis = [] #the column number of basic variables that are not xbar, x 
    row_rbasis = [] # the row number of basic slack variables (does not consider the linking rows)
    for i in 1:length(cbasis)
      if i in xbar_indices || i in x_indices
        continue
      end
      if cbasis[i] == :Basic || cbasis[i] == :NonbasicAtUpper
        push!(col_cbasis, i)
      end
    end

    row = 1
    for i in 1:length(rbasis)
      if i in link_constr_rows
        continue
      end
      if rbasis[i] == :Basic
        push!(row_rbasis, row)
      end
      row += 1
    end


    i = 1
    for col in col_cbasis
      j  = 1
      for row in 1:(constr_matrix.m)
        if row in link_constr_rows
          continue
        end
        B[j, i] = constr_matrix[row, col]
        j += 1
      end
      B[num_struct_constr + map_col_y[col] ,i] = 1
      i += 1
    end


    for row in row_rbasis
      B[row, i] = 1
      i += 1
    end

    for j in 1:length(model.colVal)
      if j in x_indices || j in xbar_indices
        continue
      end
      if cbasis[j] == :Basic || cbasis[j] ==:NonbasicAtLower
        B[num_struct_constr+map_col_y[j],i] = 1
      else
        continue
      end
      i+= 1
    end


    i = 1
    for row in 1:(constr_matrix.m) 
      if row in link_constr_rows
        continue 
      end

      new_rhs[i] = rhs[row]

      for col in  1:constr_matrix.n 
        if map_col_y[col] == 0
          if model.colVal[col] == model.colUpper[col] 
            B_N[i, col] = - constr_matrix[row, col]
            new_rhs[i] += B_N[i, col] * model.colUpper[col]
            continue
          end
        end
        B_N[i, col] = constr_matrix[row, col]
      end
      if Char(senses[row]) == 'L'
        B_N[i, constr_matrix.n + i] = 1 
      end 
      if Char(senses[row]) == 'G'
        B_N[i, constr_matrix.n + i] = -1 
      end

      i += 1
    end

    for j in 1:constr_matrix.n 
      if j in xbar_indices || j in x_indices
        continue
      end
      B_N[i, j] = 1 
      B_N[i, constr_matrix.n+num_struct_constr+map_col_y[j]] = 1
      new_rhs[i] = model.colUpper[j]
      i += 1
    end


    # #get fractional row 
    i = 1
    Binv_B_N = inv(B) * B_N
    Binv_rhs = inv(B) * new_rhs
    for col in col_cbasis
      if (model.colCat[col] == :Bin || model.colCat[col] == :Int) && (model.colVal[col]%1) < 1-1e-5 && (model.colVal[col]%1) > 1e-5 
        #generate GMI from row 
        f0 = get_frac(Binv_rhs[i])
        f = get_frac(Binv_B_N[i, :])
        cut_coeff = zeros(Float64, constr_matrix.n+num_struct_constr+num_y)
        cut_rhs = 1.0
        for j in 1:length(f)
          if abs(f[j]) < 1e-10
            continue
          end
          if j <= constr_matrix.n && (model.colCat[j] == :Int || model.colCat[j] == :Bin )
            if f[j] <= f0
              cut_coeff[j] = f[j] / f0
            else        
              cut_coeff[j] = (1-f[j]) / (1-f0)
            end
          else
            if Binv_B_N[i, j] >0
              cut_coeff[j] = Binv_B_N[i,j]/f0 
            else
              cut_coeff[j] = -Binv_B_N[i,j]/(1-f0)
            end
          end
        end

        #project slacks back into structural variables 
        for j in 1:num_struct_constr
          slack_coeff = cut_coeff[constr_matrix.n+j]
          cut_coeff -= B_N[j, constr_matrix.n + j] * slack_coeff *B_N[j, :] 
          cut_rhs -= slack_coeff * new_rhs[j] * B_N[j, constr_matrix.n + j]
          cut_coeff[constr_matrix.n+j] = 0
        end

        j = 1
        for col in 1:constr_matrix.n
          if col in x_indices || col in xbar_indices
            continue
          end
          slack_coeff = cut_coeff[constr_matrix.n+num_struct_constr+j]
          cut_coeff[col] -= slack_coeff  
          cut_rhs -= model.colUpper[col] * slack_coeff
          cut_coeff[constr_matrix.n+num_struct_constr+j] = 0
          j+=1
        end
        #change the cut back into x 
        for col in 1:constr_matrix.n 
          if col in x_indices && model.colVal[col] == model.colUpper[col]
            cut_rhs -= cut_coeff[col] * model.colUpper[col]
            cut_coeff[col] = - cut_coeff[col]
          end
        end
        cut_coeff = cut_coeff[1:constr_matrix.n]
        expr = 0.0 
        for col in 1:constr_matrix.n 
          expr += Variable(model, col) * cut_coeff[col]
        end
        # #check if the cuts cut off the fractional solution 
        # temp_lhs = 0.0 
        # for col in 1:constr_matrix.n 
        #   temp_lhs += cut_coeff[col] * model.colVal[col]
        # end
        # println(model.colVal)
        # println(cut_coeff)
        # println(cut_rhs)
        # println("new cut ")
        # println(temp_lhs-cut_rhs)
        # println(model)
        # println(expr)
        # error("GMI is generated")
        @constraint(model, expr >= cut_rhs)
      end
      i += 1

    end
    # println("========GMI==================generated")
    status = solve(model, relaxation = true)
    temp2_obj = getobjectivevalue(model)
    # println(temp2_obj-temp_obj)
    @assert status == :Optimal
    break 
  end

  @assert status == :Optimal

  dualconstraints = node.attributes[:linkconstraints]
  λnode = getdual(dualconstraints)
  nodebound = getobjectivevalue(model)

  node.attributes[:bound] = nodebound
  node.attributes[:λ] = λnode

  return status
end


function solvenodemodel(node::PlasmoNode,graph::PlasmoGraph)
  if graph.attributes[:is_nonconvex] && in_degree(graph, node) != 0
    model = node.attributes[:ubsub]
    # println(model)
    status = solve(model)

    node.attributes[:preobjval] = getvalue(node.attributes[:preobj])
  else 
    model = getmodel(node)
    status = solve(model)
    if status != :Optimal && in_degree(graph, node) == 0
      graph.attributes[:is_infeasible] = true
      graph.attributes[:LB] = +Inf
    end
    if status != :Optimal && in_degree(graph, node) !=0
      println(node.attributes[:xin])
      error("upper bound subproblem not solved to optimality")
    end
    if in_degree(graph,node) == 0 # Root node
      graph.attributes[:LB] = getobjectivevalue(model)
    end
    if graph.attributes[:LB] > graph.attributes[:global_UB]
      graph.attributes[:fathomed_by_bound] = true
    end
    node.attributes[:preobjval] = getvalue(node.attributes[:preobj])
  end

end

function takex(node::PlasmoNode, graph::PlasmoGraph)
  if in_degree(graph, node) == 0 
    return true
  end
  xinvals = node.attributes[:xin]
  # println(xinvals)
  xinvars = node.attributes[:xinvars]
  if length(xinvals) > 0
    fix.(xinvars,xinvals)
    if graph.attributes[:is_nonconvex]
      fix.(node.attributes[:ubxinvars], xinvals)
    end
  end
end

function putx(node::PlasmoNode,graph::PlasmoGraph)
  childvars = node.attributes[:childvars]
  children = out_neighbors(graph,node)
  length(children) == 0 && return true

  for child in children
    xnode = getvalue(childvars[getnodeindex(graph,child)])
    child.attributes[:xin] = xnode
  end
end

function putcutdata(node::PlasmoNode,graph::PlasmoGraph,cuts::Array{Symbol,1})
  parents = in_neighbors(graph,node)
  length(parents) == 0 && return true
  parent = parents[1]    # Assume only one parent
  parentcuts = parent.attributes[:cutdata]
  θk = node.attributes[:bound]
  λk = node.attributes[:λ]
  xk = node.attributes[:xin]
  nodeindex = getnodeindex(graph,node)
  if :LP in cuts || :GMI in cuts
    bcd = BendersCutData(θk, λk, xk)
    push!(parentcuts[nodeindex],bcd)
  end
  if :LLinteger in cuts
    llintcd = LLIntegerCutData(θlb,xk)
    push!(parentcuts[nodeindex],llintcd)
  end
  if :Integer in cuts
    intcd = IntegerCutData(xk)
    push!(parentcuts[nodeindex],intcd)
  end
end

function generatecuts(node::PlasmoNode,graph::PlasmoGraph)
  children = out_neighbors(graph,node)
  length(children) == 0 && return true

  cutdataarray = node.attributes[:cutdata]
  previouscuts = node.attributes[:prevcuts]
  thisitercuts = Dict()
  samecuts = Dict()
  for child in children
    childindex = getnodeindex(graph,child)
    thisitercuts[childindex] = CutData[]
    samecuts[childindex] = Bool[]

    while length(cutdataarray[childindex]) > 0
      cutdata = pop!(cutdataarray[childindex])
      samecut = in(cutdata,previouscuts[childindex])
      push!(samecuts[childindex],samecut)
      samecut && continue
      if typeof(cutdata) == BendersCutData
        generatebenderscut(node,cutdata,childindex)
      elseif typeof(cutdata) == LLIntegerCutData
        generateLLintegercut(node,cutdata)
      elseif typeof(cutdata) == IntegerCutData
        generateintegercut(node,cutdata)
      end
      push!(thisitercuts[childindex],cutdata)
    end
    samecuts[childindex] = reduce(*,samecuts[childindex]) && length(samecuts[childindex]) > 0
  end
  node.attributes[:prevcuts] = thisitercuts
  nodesamecuts = collect(values(samecuts))
  node.attributes[:stalled] = reduce(*,nodesamecuts)
  node.attributes[:stalled] && warn("Node $(node.label) stalled")
  if in(node,graph.attributes[:roots]) && node.attributes[:stalled]
    graph.attributes[:stalled] = true
 end
end

function generatebenderscut(node::PlasmoNode, cd::BendersCutData,index)
  model = getmodel(node)
  θ = getindex(model, :θ)
  x = node.attributes[:childvars][index]
  @constraint(model, θ[index] >= cd.θk + cd.λk'*(cd.xk - x))
end


function identifylevels(graph::Plasmo.PlasmoGraph)
  #Create lists of root and leaf nodes in graph
  roots = graph.attributes[:roots] = []
  leaves = graph.attributes[:leaves] = []
  #Create dictionary to keep track of levels of nodes
  levels = graph.attributes[:levels] = Dict()
  #Iterate through every node to check for root/leaf nodes
  for node in values(graph.nodes)
    node.attributes[:xin] = []
    node.attributes[:λ] = []
    node.attributes[:bound] = NaN
    node.attributes[:xinvars] = []
    node.attributes[:preobjval] = NaN
    node.attributes[:linkconstraints] = []
    if graph.attributes[:is_nonconvex] && in_degree(graph, node) !=0 
      node.attributes[:ubxinvars] = []
    end
    #If the node does not have parents it is a root node
    if in_degree(graph,node) == 0
      push!(roots,node)
    end
    #If the node does not have children it is a leaf node
    if out_degree(graph,node) == 0
      push!(leaves,node)
    end
  end

  #Start mapping level from the root nodes
  current = roots
  level = 1
  levels[1]  = roots
  while length(current) > 0
    levels[level] = current
    children = []
    for node in current
        push!(children,out_neighbors(graph,node)...)
        node.attributes[:childvars] = Dict(getnodeindex(graph,child) => [] for child in out_neighbors(graph,node))
        node.attributes[:cutdata] = Dict(getnodeindex(graph,child) => CutData[] for child in out_neighbors(graph,node))
        node.attributes[:prevcuts] = Dict(getnodeindex(graph,child) => CutData[] for child in out_neighbors(graph,node))
        node.attributes[:stalled] = false
    end
    current = children
    level += 1
  end
  graph.attributes[:numlevels] = level - 1
end

function bdprepare(graph::Plasmo.PlasmoGraph, cuts::Array{Symbol,1}=[:LP])
  if haskey(graph.attributes,:preprocessed)
    return true
  end
  graph.attributes[:solution]= Solution(method=:benders)
  identifylevels(graph)
  graph.attributes[:normalized] = normalizegraph(graph)
  graph.attributes[:stalled] = false
  graph.attributes[:mflat] = create_flat_graph_model(graph)
  graph.attributes[:UB] = Inf
  graph.attributes[:global_UB] = +Inf
  graph.attributes[:fathomed_by_bound] = false
  graph.attributes[:is_infeasible] = false
  setsolver(graph.attributes[:mflat],graph.solver)

  links = getlinkconstraints(graph)
  numlinks = length(links)

  for index in 1:length(graph.nodes)
    node = graph.nodes[index]
    model = getmodel(node)
    if model.solver == JuMP.UnsetSolver()
      model.solver = graph.solver
    end
    if graph.attributes[:is_nonconvex] && in_degree(graph, node) != 0
      node.attributes[:preobj] = node.attributes[:ubsub].obj 
    else
      node.attributes[:preobj] = model.obj
    end

    #create GMI if :GMI is in cuts
    if :GMI in cuts
      for col in 1:length(model.colUpper)
        if model.colUpper[col] > 1e10 
          setupperbound(Variable(model, col), 1e10)
        end
      end
    	node.attributes[:GMI] = deepcopy(model.internalModel)
      println("copy")
    end

    #map variable name to variable 
    node.attributes[:varname_to_var] = Dict()
    node.attributes[:varname_to_varcat] = Dict()
    for i in 1:length(model.colNames)
      varname = model.colNames[i]
      node.attributes[:varname_to_var][varname] = Variable(model, i)
      node.attributes[:varname_to_varcat][varname] = model.colCat[i]
    end

    if graph.attributes[:is_nonconvex] && in_degree(graph, node) != 0
      node.attributes[:nonconvex_varname_to_var] = Dict()
      for i in 1:length(node.attributes[:ubsub].colNames)
        varname = node.attributes[:ubsub].colNames[i]
        node.attributes[:nonconvex_varname_to_var][varname] = Variable(node.attributes[:ubsub], i)
      end
    end

    if out_degree(graph,node) != 0
#Add theta to parent nodes
      childrenindices = [getnodeindex(graph,child) for child in out_neighbors(graph,node)]
      sort!(childrenindices)
      @variable(model, θ[i in childrenindices] >= -1e6)
      model.obj += sum(θ[i] for i in childrenindices)
      node.attributes[:θ_to_scenario] = Dict()
      node.attributes[:scenario_to_θ] = Dict()
      for child in out_neighbors(graph,node)
        if haskey(child.attributes, :scenario)
          index = getnodeindex(graph, child)
          node.attributes[:θ_to_scenario][index] = child.attributes[:scenario]
          node.attributes[:scenario_to_θ][child.attributes[:scenario]] = index
        end
      end

    end
  end
  #Add dual constraint to child nodes using the linking constraints
  for (numlink,link) in enumerate(links)
    #Take the two variables of the constraint
    var1 = link.terms.vars[1]
    var2 = link.terms.vars[2]
    #Determine which nodes they belong to
    nodeV1 = getnode(var1)
    nodeV2 = getnode(var2)
    #Set the order of the nodes
    if ischildnode(graph,nodeV1,nodeV2)
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
    childindex = getnodeindex(graph,childnode)
    childmodel = getmodel(childnode)
    push!(parentnode.attributes[:childvars][childindex],parentvar)
    valbar = @variable(childmodel)
    setname(valbar,"varlbar$numlink")
    push!(childnode.attributes[:xinvars],valbar)
    conref = @constraint(childmodel, valbar - childvar == 0)
    push!(childnode.attributes[:linkconstraints], conref)
  end

  #link master and ub subproblem for nonconvex 
  if graph.attributes[:is_nonconvex]
    for (numlink,link) in enumerate(links)
      #Take the two variables of the constraint
      var1 = link.terms.vars[1]
      var2 = link.terms.vars[2]
      #Determine which nodes they belong to
      nodeV1 = getnode(var1)
      nodeV2 = getnode(var2)
      #Set the order of the nodes
      if ischildnode(graph,nodeV1,nodeV2)
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
      childindex = getnodeindex(graph,childnode)
      childmodel = childnode.attributes[:ubsub]
      # println(childmodel)
      # println(childmodel.colNames)
      temp_model = childvar.m 
      childvar_name = temp_model.colNames[childvar.col]
      # println(childvar_name)
      # println(get_col_from_varname(childmodel, childvar_name))
      # push!(parentnode.attributes[:childvars][childindex],parentvar)
      valbar = @variable(childmodel)
      setname(valbar,"varlbar$numlink")
      push!(childnode.attributes[:ubxinvars],valbar)
      @constraint(childmodel, valbar - Variable(childmodel,get_col_from_varname(childmodel, childvar_name)) == 0)
    end
  end
  
  
  graph.attributes[:preprocessed] = true
end


