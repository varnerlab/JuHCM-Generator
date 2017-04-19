# my files -
include("DataDictionary.jl")
include("Balances.jl")
include("Control.jl")
include("Kinetics.jl")
include("Dilution.jl")
include("SolveBalances.jl")

# we are using SUNDIALS to solve the balances -
using Sundials
using PyPlot
