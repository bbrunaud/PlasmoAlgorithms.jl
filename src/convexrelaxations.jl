function add_McCormick(m::Model, x::Variable, y::Variable, xy::Variable)
	x_ub = getupperbound(x)
	x_lb = getlowerbound(x)
	y_ub = getupperbound(y)
	y_lb = getlowerbound(y)
	@constraint(m, xy >= x_lb * y + x * y_lb - x_lb * y_lb)
	@constraint(m, xy>= x_ub * y + x * y_ub - x_ub * y_ub)
	@constraint(m, xy<= x_ub * y + x * y_lb - x_ub * y_lb)
	@constraint(m, xy<= x_lb * y + y_ub * x - x_lb * y_ub)
end

# the x variable will be discretize over its domain. n is the number of intervals
function add_PiecewiseMcCormick(m::Model, x::Variable, y::Variable, xy::Variable, n::Int64, xmap::Dict)
	x_ub = getupperbound(x)
	x_lb = getlowerbound(x)
	y_ub = getupperbound(y)
	y_lb = getlowerbound(y)
	x_points = ones(n+1) * x_lb
	for i in 1:n
		x_points[i+1] = x_lb + i * (x_ub-x_lb) / n 
	end
	#find the binary variables for piecewise linear representation 
	delta = []
	xname = getname(x)
	yname = getname(y)
	if haskey(xmap, x)
		delta = xmap[x][:delta]
		dot_x = xmap[x][:dot_x]
	else
		delta = @variable(m, delta[i in 1:n], Bin, basename="delta_$xname")
		xmap[x] = Dict()
		xmap[x][:delta] = delta
		@constraint(m, sum(delta[i] for i in 1:n) == 1)
		dot_x = @variable(m, dot_x[i in 1:n], basename="dot_$xname")
		xmap[x][:dot_x] = dot_x
		@constraint(m, sum(dot_x[i] for i in 1:n) == x)
		@constraint(m, [i in 1:n], dot_x[i]>= x_points[i] * delta[i])
		@constraint(m, [i in 1:n], dot_x[i] <= x_points[i+1] * delta[i])
	end
	@variable(m, dot_xy[i in 1:n], basename="dot$xname$yname")
	@variable(m, dot_y[i in 1:n], basename="dot$yname")
	@constraint(m, sum(dot_xy[i] for i in 1:n) == xy)
	@constraint(m, sum(dot_y[i] for i in 1:n) == y)
	@constraint(m, [i in 1:n], dot_y[i] <= y_ub * delta[i])
	@constraint(m, [i in 1:n], dot_y[i] >= y_lb * delta[i])
	@constraint(m, [i in 1:n], dot_xy[i] >= x_points[i] * dot_y[i] + dot_x[i] * y_lb - x_points[i] * y_lb * delta[i])
	@constraint(m, [i in 1:n], dot_xy[i]>= x_points[i+1] * dot_y[i] + dot_x[i] * y_ub - x_points[i+1] * y_ub * delta[i])
	@constraint(m, [i in 1:n], dot_xy[i]<= x_points[i+1] * dot_y[i] + dot_x[i] * y_lb - x_points[i+1] * y_lb * delta[i])
	@constraint(m, [i in 1:n], dot_xy[i]<= x_points[i] * dot_y[i] + y_ub * dot_x[i] - x_points[i] * y_ub * delta[i])
end




# the x variable will be discretize over its domain. n is the number of binary digits
#xmap is used to record whether x has been discretized before 
#the implementation corresponds to equation 16 in Misener and Floudas (2011)
#APOGEE: Global optimization of standard, generalized, and extended pooling problems via linear and logarithmic partitioning schemes
function add_LogPiecewiseMcCormick(m::Model, x::Variable, y::Variable, xy::Variable, n::Int64, xmap::Dict)
	x_ub = getupperbound(x)
	x_lb = getlowerbound(x)
	y_ub = getupperbound(y)
	y_lb = getlowerbound(y)
	a = (x_ub - x_lb) / (2^n)
	#find the binary variables for logrithmic piecewise  representation 
	lambda = []
	xname = getname(x)
	yname = getname(y)
	if haskey(xmap, x)
		lambda = xmap[x][:lambda]
	else
		lambda = @variable(m, lambda[i in 1:n], Bin,  basename="lambda_$xname")
		@constraint(m, x_lb +sum(2^(i-1) * a * lambda[i] for i in 1:n) <= x )
		@constraint(m, x_lb +sum(2^(i-1) * a * lambda[i] for i in 1:n) + a >= x )
		xmap[x] = Dict()
		xmap[x][:lambda] = lambda
	end
	@variable(m, s[i in 1:n]>=0, basename="s_$xname$yname")
	@variable(m, dot_y[i in 1:n]>=0, basename="dot$yname$xname")

	@constraint(m, [i in 1:n], dot_y[i] <= (y_ub - y_lb) * lambda[i])
	@constraint(m, [i in 1:n], dot_y[i] == (y - y_lb) - s[i] )
	@constraint(m, [i in 1:n], s[i] <= (y_ub - y_lb) * (1 - lambda[i]) )
	@constraint(m, xy >= x * y_lb + x_lb * (y - y_lb) + sum(a * 2^(i-1) * dot_y[i] for i in 1:n))
	@constraint(m, xy >= x*y_ub + (x_lb+a)*(y- y_ub) + sum(a* 2^(i-1) * (dot_y[i] - (y_ub - y_lb) * lambda[i]) for i in 1:n))
	@constraint(m, xy <= x * y_lb + (x_lb + a) * (y - y_lb) + sum(a * 2^(i-1) * dot_y[i] for i in 1:n))
	@constraint(m, xy <= x*y_ub + x_lb*(y- y_ub) + sum(a* 2^(i-1) * (dot_y[i] - (y_ub - y_lb) * lambda[i]) for i in 1:n))
end



#This is convex relaxation for bilinear programs. We construct the hull relaxation of  CNF for piecewise McCormick
function add_disjunction(model::Model, A_lp,b_lp,m_lp, n, x_var, y_var_array, xy_var_array, npoints)
    lambda = @variable(model, 0<=lambda[i in 1:npoints]<=1, basename="lambda_$x_var")
    @constraint(model, sum(lambda[i] for i in 1:npoints) == 1)
    x_ub = getupperbound(x_var)
    x_lb = getlowerbound(x_var)
    interval_length = (x_ub - x_lb) / npoints
    x_points = ones(npoints+1) * x_lb
    for d in 1:npoints
        x_points[d+1] = x_lb + interval_length * d 
    end 
    x_index = x_var.col 
    y_index = []
    xy_index = []
    for i in 1:length(y_var_array)
        push!(y_index, y_var_array[i].col)
        push!(xy_index, xy_var_array[i].col)
    end

    # #add constraints for each disjunction 
    vars_all_djc = []
    for d in 1:npoints
        var_this_djc =[] 
        for i in 1:n
            temp = @variable(model, basename="dot_" * model.colNames[i] * "_$x_var$d")
            push!(var_this_djc, temp)
        end

        @constraint(model, [i in 1:m_lp], sum(A_lp[i, j] * var_this_djc[j] for j in 1:n) >= b_lp[i] * lambda[d])
        #add McCormick for vars in y_var_array
        dot_x = var_this_djc[x_index]
        for i in 1:length(y_var_array)
            dot_y = var_this_djc[y_index[i]]
            y = y_var_array[i]
            dot_xy = var_this_djc[xy_index[i]]
            xy = xy_var_array[i]
            y_ub = getupperbound(y)
            y_lb = getlowerbound(y)
            @constraint(model,  dot_xy >= x_points[d] * dot_y + dot_x * y_lb - x_points[d] * y_lb * lambda[d])
            @constraint(model,  dot_xy>= x_points[d+1] * dot_y + dot_x * y_ub - x_points[d+1] * y_ub * lambda[d])
            @constraint(model, dot_xy<= x_points[d+1] * dot_y + dot_x * y_lb - x_points[d+1] * y_lb * lambda[d])
            @constraint(model, dot_xy<= x_points[d] * dot_y + y_ub * dot_x - x_points[d] * y_ub * lambda[d])
        end

        push!(vars_all_djc, var_this_djc)
    end
    @constraint(model, [i in 1:n], Variable(model, i) == sum(vars_all_djc[d][i] for d in 1:npoints))
end

function getModelInfo(model::Model)
	A = JuMP.prepConstrMatrix(model)
    n = size(A,2)                               # Numver of variables of original problem
    m = size(A,1)                               # Number of constraints of original problem
    conBounds = JuMP.prepConstrBounds(model)       # trabsforn all constraints in greater than
    LHS = conBounds[1]
    RHS = conBounds[2]
    b = zeros(m)
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

    vars = Variable.(model, 1:n)            # Get a vector of all the variables

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

    return A_lp,b_lp,m_lp, n

end












