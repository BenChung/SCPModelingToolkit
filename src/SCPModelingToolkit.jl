module SCPModelingToolkit
using ModelingToolkit, SCPToolbox
import Symbolics


include("problem.jl")
include("utilities.jl")
include("variables.jl")
include("adapter.jl")
include("initialization.jl")

end
