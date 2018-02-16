mutable struct Solution
    problemname
    method::Symbol
    solvetime
    objval::Float64
    bestbound::Float64
    gap::Float64
    numiterations::Int64
    termination
    # Iteration Data
    iterval::Array{Float64,1}
    iterbound::Array{Float64,1}
    itertime::Array{Float64,1} # in seconds
    # Lagrange
    α
    step
end

function Solution(;method=:none)
    return Solution("",method,0,NaN,NaN,NaN,0,"Not Executed",Float64[],Float64[],Float64[],Float64[],Float64[])
end

function saveiteration(s::Solution,tstamp::Float64,arr::Array{Float64,1},n=1)
    push!(s.iterval,arr[1])
    push!(s.iterbound,arr[2])
    push!(s.itertime,arr[3])
    if length(arr) == 5 && s.method  == :dual_decomposition
        push!(s.α,arr[4])
        push!(s.step,arr[5])
    end
    if n == 1
        s.bestbound = maximum(s.iterbound)
        s.objval = minimum(s.iterval)
    else
        s.bestbound = minimum(s.iterbound)
        s.objval = maximum(s.iterval)
    end
    s.numiterations += 1
    s.gap = abs(s.objval - s.bestbound)/s.objval
    s.solvetime = tstamp
end

function Base.show(io::IO,s::Solution)
    println(" SOLUTION SUMMARY")
    println("-----------------")
    println("Method: $(s.method)")
    println("Objective Value : $(s.objval)")
    println("Best Bound : $(s.bestbound)")
    println("Gap : $(round(s.gap*100,2)) %")
    println("Iterations : $(s.numiterations)")
    println("Solution Time : $(round(s.solvetime,2)) s")
    println("Termination : $(s.termination)")
end
