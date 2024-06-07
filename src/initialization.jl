
struct StraightlineInterpolate <: Initialization
    initial::Dict{SymbolicUtils.BasicSymbolic, Float64}
    final::Dict{SymbolicUtils.BasicSymbolic, Float64}
    parameters::Dict{SymbolicUtils.BasicSymbolic, Float64}
    StraightlineInterpolate(
        initial::Dict{SymbolicUtils.BasicSymbolic, Float64},
        final::Dict{SymbolicUtils.BasicSymbolic, Float64},
        parameters::Dict{SymbolicUtils.BasicSymbolic, Float64}) = new(initial, final, parameters)
    StraightlineInterpolate(
        initial::Dict{Num, Float64},
        final::Dict{Num, Float64},
        parameters::Dict{Num, Float64}) = 
        new(unwrap_dict(initial), unwrap_dict(final), unwrap_dict(parameters))
end

function initialize(pbm::TrajectoryProblem, p::SCPtProblem, s::StraightlineInterpolate)
    ikeys = keys(s.initial)
    fkeys = keys(s.final)
    akeys = union(intersect(ikeys, fkeys), keys(s.parameters))
    @assert all(s -> s ∈ akeys, p.states) && all(s -> s ∈ akeys, p.controls) && all(s -> s ∈ akeys, p.parameters) "All states, controls, and parameters must be have initial and final values provided"

    x0 = ModelingToolkit.varmap_to_vars(s.initial, p.states)
    u0 = ModelingToolkit.varmap_to_vars(s.initial, p.controls)
    pv = ModelingToolkit.varmap_to_vars(s.parameters, p.parameters)
    xf = ModelingToolkit.varmap_to_vars(s.final, p.states)
    uf = ModelingToolkit.varmap_to_vars(s.final, p.controls)
    
    problem_set_guess!(pbm, (N, pbm) -> begin
        res = straightline_interpolate(x0, xf, N), straightline_interpolate(u0, uf, N), isnothing(pv) ? Float64[] : pv
        return res
    end)
end