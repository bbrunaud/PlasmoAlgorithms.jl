function solvegmirelaxation(node::PlasmoNode)
  # model = node.attributes[:GMI]
  model = getmodel(node)
  status=solve(model, relaxation=true)
  temp_obj = getobjectivevalue(model)
  while true 
    cpx = model.internalModel.inner 


    constr_matrix = CPLEX.get_constr_matrix(cpx)
    rhs = CPLEX.get_rhs(cpx)
    #identify xbar x, linking constraints
    link_constr_rows = []
    xbar_indices = []
    xbar_value = []
    xbar_ub = []
    xbar_lb = []
    x_indices = []

    num_linkconstr = length(node.attributes[:linkconstraints])
    idx = node.attributes[:linkconstraints][1].idx 
    for i in 1:num_linkconstr
      row = constr_matrix[i+idx-1,:]
      # is_one = false 
      # is_negative_one = false 
      # if length(nonzeros(row)) == 2
      #   for j in 1:2
      #     if nonzeros(row)[j] == 1
      #       is_one = true 
      #     end
      #     if nonzeros(row)[j] == -1 
      #       is_negative_one = true
      #     end
      #   end
      # end
      push!(link_constr_rows, i + idx-1)
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

    map_col_y = zeros(Int64, length(model.colVal)) # give a column return it's the i th y, y represent the second stage decisions
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

    #get the basic and nonbasic variables
    cbasis = CPLEX.get_basis(cpx)[1]
    rbasis = CPLEX.get_basis(cpx)[2]
    senses = CPLEX.get_constr_senses(cpx)

    num_struct_constr = constr_matrix.m - num_linkconstr
    num_y = length(model.colVal) - 2*num_linkconstr
    num_x = num_linkconstr
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
    println("========GMI==================generated")
    status = solve(model, relaxation = true)
    temp2_obj = getobjectivevalue(model)
    println(temp2_obj-temp_obj)
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