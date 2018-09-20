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




















