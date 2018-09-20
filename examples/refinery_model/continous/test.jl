using JuMP
using PlasmoAlgorithms
using Plasmo
include("input.jl")
include("lagsub.jl")
include("benderssub.jl")


model = generate_benderssub(prob=prob[1], Crude_yield_data = Crude_yield_data[:,:,1], Desulphurisation_cost=Desulphurisation_cost[:,1], Sulphur_2=Sulphur_2[:,1], Sulphur_GO_data= Sulphur_GO_data[:,1])
































