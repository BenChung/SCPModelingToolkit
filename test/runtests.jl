using SCPModelingToolkit, LinearAlgebra, ModelingToolkit
using SCPToolbox
import Symbolics, ECOS

using Test

@testset "SCPModelingToolkit.jl" begin
    import("scpt_examples/quad.jl")
end
