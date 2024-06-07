module SCPModelingToolkit
using ModelingToolkit, SCPToolbox, ForwardDiff
import Symbolics, JuMP


include("problem.jl")
include("utilities.jl")
include("variables.jl")
include("adapter.jl")
include("initialization.jl")

end
