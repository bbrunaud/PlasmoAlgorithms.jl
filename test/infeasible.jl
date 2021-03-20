using JuMP
using Gurobi
using Xpress

optimizer = Gurobi.Optimizer
#optimizer = Xpress.Optimizer

m = Model()
JuMP.set_optimizer(m, optimizer)
set_optimizer_attribute(m, "InfUnbdInfo", 1)
@variable(m, x >= 0)
@variable(m, y >= 0)
@constraint(m, c1, x + y <= 10)
@constraint(m, c2, x == 12)
@objective(m, Max, x + y)

optimize!(m)

print(dual_status(m))
println("Dual c1: $(dual(c1))")
println("Dual c1: $(dual(c2))")
