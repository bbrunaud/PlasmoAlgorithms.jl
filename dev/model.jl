using CPLEX
using JuMP
function generate_model()
	m = Model(solver=CplexSolver())
	@variable(m, 0<=x1<=1)
	@variable(m, 0<=x2<=1)
	@variable(m, 1<=x1bar<=1)
	@variable(m, 1<=x2bar<=1)
	@variable(m, y1, Bin)
	@variable(m, y2, Bin)
	@variable(m, y3, Bin)
	@variable(m, y4, Bin)
	@variable(m, 0<=R<=1000, Int)

	@constraint(m, x1==x1bar)
	@constraint(m, x2==x2bar)
	@constraint(m, x1+2*y1+3*y2+4*y3+5*y4-R<= 10)
	@constraint(m, x2+6*y1+y2+3*y3+2*y4-R<=3)
	@objective(m, Min, -(16*y1+19*y2+23*y3+28*y4-100*R))
	return m 
end