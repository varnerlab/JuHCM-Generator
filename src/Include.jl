# required external packages -
using ArgParse
using JSON

# components of JuHCM -
include("Types.jl")
include("Macros.jl")
include("Parser.jl")
include("Problem.jl")
include("./strategy/JuliaStrategy.jl")
include("Common.jl")
