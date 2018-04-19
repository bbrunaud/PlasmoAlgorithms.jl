include("ft2level.jl")

iter = 1
bval = 200

ms2 = deepcopy(ms)

for i in 1:5
run(`cowsay Iteration $iter - b = $bval`)

m = opmodel(bval)
solve(m, relaxation= true)
lrd = dual(m)
lprel = getobjectivevalue(m)

@constraint(ms, θ >= lprel - lrd*(bval - b))
solve(ms)

bval = getvalue(b)
iter += 1
end

iter = 1
bval = 200

for i in 1:4
run(`cowsay Iteration $iter - b = $bval`)

m = opmodel(bval)
n = opmodel(bval,:NLP)
solve(m)
n.colVal = m.colVal
solve(n)
nd = dual(n)
nobj = getobjectivevalue(n)

θ = getindex(ms2,:θ)
b = getindex(ms2,:b)
@constraint(ms2, θ >= nobj - nd*(bval - b))



solve(ms2)

bval = getvalue(b)
iter += 1
end
