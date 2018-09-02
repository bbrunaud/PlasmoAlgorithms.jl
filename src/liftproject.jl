function getmatrixform(node::PlasmoNode)
    # get info on the original model
    sp = getmodel(node)
    A = JuMP.prepConstrMatrix(sp)
    n = size(A,2)                               # Numver of variables of original problem
    m = size(A,1)                               # Number of constraints of original problem
    conBounds = JuMP.prepConstrBounds(sp)       # trabsforn all constraints in greater than
    LHS = conBounds[1]
    RHS = conBounds[2]
    b = zeros(m)


    # Re-arrenge Matrix in the canonical form Ax>=b
    # TO DO : handle equality constraits
    eq_counter = 0
    for i = 1:m
        if LHS[i] == RHS[i]         # This identifies an equality constraint
            eq_counter = eq_counter +1

            # Note that the order of the assignement here is important
            # otherwise you would change the sign of the coefficient of the A matrix
            A = [A;A[i,:]']
            b = [b;LHS[i]]

            A[i,:] = -A[i,:]        # Add one inequality
            b[i] = -LHS[i]



        else                        # This is inequality constraint
            if RHS[i] !=Inf         # This identifies a less or equal constraint
                b[i] = -RHS[i]
                A[i,:] = -A[i,:]
            else
                b[i] = LHS[i]
            end
        end
    end

    # Update matrix formulation
    if eq_counter != 0
        m = size(A,1)              # Update Number of inequalities
    end
    node.attributes[:A] = A
    node.attributes[:b] = b
    node.attributes[:m] = m
    node.attributes[:n] = n 

end

#create CGLPs for each child node. Need to call this function again whenever branching is applied on this node 
function setCGLP(node::PlasmoNode, graph::PlasmoGraph)

    #retrieve the original matrix 
    A = node.attributes[:A]
    b = node.attributes[:b]
    m = node.attributes[:m]
    n = node.attributes[:n]

    sp = getmodel(node)
    # Create a map from varaibles to binary variables
    vars = Variable.(sp, 1:n)            # Get a vector of all the variables

    # This dictionary contains the binary variables and their position (columsn of A)
    node.attributes[:map_binary] = Dict()
    for i = 1:n
        varCategory = getcategory(vars[i])
        if  varCategory == :Bin && (!(i in node.attributes[:linking_vars_indices]))
            node.attributes[:map_binary][vars[i]] = i
        elseif varCategory == :Int
            println("==================================================")
            println("WARNING: The subproblems present Integer Variables")
            println("WARNING: The Lift-and-Project cuts will not define a facet inequality, therefore the convex hull is not correctly identified")
            println("WARNING: It is very likely the solution will not be optimal")
            println("==================================================")
        end
    end
    map_binary = node.attributes[:map_binary]

    # binary Info
    num_binary = length(map_binary)

    #create LP relaxation matrix including bounds of all the variables 
    m_lp = m 
    b_lp = b 
    A_lp = A
    for i in 1:n 
        ub = getupperbound(vars[i])
        lb = getlowerbound(vars[i])
        if ub < 1e6
            new_row = zeros(n)'
            new_row[i] = -1 
            A_lp = [A_lp; new_row]
            b_lp = [b_lp; -ub]
            m_lp += 1
        end

        if lb > -1e6
            new_row = zeros(n)'
            new_row[i] = 1
            A_lp = [A_lp; new_row]
            b_lp = [b_lp; lb]
            m_lp += 1
        end
    end

    #create CGLPs 
    node.attributes[:CGLPs] =  Dict()
    for var in keys(map_binary)
        CGLP = Model()
        CGLP.solver = graph.solver
        @variables CGLP begin
            α[1:n]
            β
            u[1:m_lp] >= 0
            v[1:m_lp] >= 0
            u0 >= 0
            v0 >= 0
        end     

        # CGLP formulation (Balas)
        ek = zeros(n)
        ek[map_binary[var]] = 1         # Pick the disjunction on the j fractional variable
        # Convex Hull reformulation
        @constraint(CGLP, convexDual1[row in 1:n], α[row] - sum(u[ii]*A_lp[ii,row]  for ii in 1:m_lp) + u0*ek[row] == 0)
        @constraint(CGLP, convexDual2[row in 1:n], α[row] - sum(v[ii]*A_lp[ii,row]  for ii in 1:m_lp) - v0*ek[row] == 0)
        @constraint(CGLP, convexDual3, -β + sum(u[ii]*b_lp[ii] for ii in 1:m_lp)      == 0)
        @constraint(CGLP, convexDual4, -β + sum(v[ii]*b_lp[ii] for ii in 1:m_lp) + v0 == 0)
        # normalization constraint
        @constraint(CGLP, normalization, sum(u[ii] + v[ii] for ii in 1:m_lp) + u0 + v0 == 1)   
        node.attributes[:CGLPs][map_binary[var]] = CGLP
    end

end


function solveliftandprojectrelaxation(node::PlasmoNode, graph::PlasmoGraph)
    model = getmodel(node)
    status = solve(model, relaxation=true)
    #only start adding lift and project cuts when the LP relaxation has stalled 
    if graph.attributes[:roots][1].attributes[:LP_stalled]
        #get the values of variables
        n = node.attributes[:n]
        vars = Variable.(model, 1:n)
        x_bar = getvalue(vars)
        num_frac = 0 
        for index in values(node.attributes[:map_binary])
            if abs(x_bar[index]%1)>1e-3 && abs(x_bar[index]%1) < 1-1e-3
                num_frac += 1
                CGLP = node.attributes[:CGLPs][index]
                α = getindex(CGLP, :α)
                β = getindex(CGLP, :β)
                @objective(CGLP, Min, sum(α[ii]*x_bar[ii] for ii in 1:n)  - β)
                # Solve the Cut Generation problem
                CGLPstatus = solve(CGLP)
                if CGLPstatus == :Optimal
                    #add cuts to model 
                    α_cut = getvalue(α)
                    β_cut = getvalue(β)
                    @constraint(model, sum(α_cut[i]*vars[i]  for i in 1:n) >= β_cut)
                else
                    println("ERROR: ")
                    println("ERROR: Cut Generation Linear problem not Optimal")
                    println("ERROR: It is not possible to generate Lift-and-Project cuts")
                    error("lp cuts error")
                    break
                end
            end
        end

        if num_frac >0 
            solve(model, relaxation = true)
        end
    end

    dualconstraints = node.attributes[:linkconstraints]
    λnode = getdual(dualconstraints)
    nodebound = getobjectivevalue(model)

    node.attributes[:bound] = nodebound
    node.attributes[:λ] = λnode
end






















