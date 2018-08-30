using PlasmoAlgorithms
using JLD
using DataFrames


d = load("Standard.jld","d")

df = DataFrame(Item = [],
               P0 = [],
               T0 = [],
               PT0 = [],
               PR = [],
               TR = [],
               PTR = [])

res = zeros(6,6)
for (k,s) in enumerate([6 20])
    for (j,rel) in enumerate([:zero,:relaxation])
        for (i,p) in enumerate(["P", "T", "PT"])
            r = d["$p$s.jl",:subgradient,rel]
            
            res[3k-2:3k,i+3*(j-1)] = [r.numiterations, r.solvetime/r.numiterations,r.gap*100]
            println("$p $s $rel : $(3k),($i+3*($j-1))")
        end
    end
end

