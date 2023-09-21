
struct ConicConstraint
    name::String
    eqn::Vector{SymbolicUtils.BasicSymbolic}
    cone::SCPToolbox.Parser.ConicLinearProgram.SupportedCone
    ConicConstraint(name::String, eqn::SymbolicUtils.BasicSymbolic, cone::SCPToolbox.Parser.ConicLinearProgram.SupportedCone) = new(name, [eqn], cone)
    ConicConstraint(name::String, eqn::Num, cone::SCPToolbox.Parser.ConicLinearProgram.SupportedCone) = new(name, [Symbolics.unwrap(eqn)], cone)
    ConicConstraint(name::String, eqns::Vector{SymbolicUtils.BasicSymbolic}, cone::SCPToolbox.Parser.ConicLinearProgram.SupportedCone) = new(name, eqns, cone)
    ConicConstraint(name::String, eqns::Vector{Num}, cone::SCPToolbox.Parser.ConicLinearProgram.SupportedCone) = new(name, Symbolics.unwrap.(eqns), cone)
    
end
abstract type Initialization end

struct SCPtProblem 
    states::Vector{SymbolicUtils.BasicSymbolic}
    controls::Vector{SymbolicUtils.BasicSymbolic}
    parameters::Vector{SymbolicUtils.BasicSymbolic}
    scale_advice::Dict{SymbolicUtils.BasicSymbolic, Tuple{Float64, Float64}}

    dynamics::ModelingToolkit.ODESystem
    constraints::Vector{ConicConstraint}
    terminal_cost::Union{SymbolicUtils.BasicSymbolic, Float64}
    running_cost::Function
    nonconvex_constraints::Vector{SymbolicUtils.BasicSymbolic}
    bcs::Dict{Symbol, Vector{SymbolicUtils.BasicSymbolic}}
    initalizer::Initialization
    function SCPtProblem(; 
        scale_advice::Dict{Num, Tuple{Float64, Float64}} = Dict(), 
        dynamics::ModelingToolkit.ODESystem,
        constraints::Vector{ConicConstraint},
        terminal_cost::Union{Float64, Num},
        running_cost::Function, 
        nonconvex_constraints::Vector{Num},
        bcs::Dict{Symbol, Vector{Num}},
        initalizer::Initialization)
        vars = unique([ModelingToolkit.states(dynamics);
                ModelingToolkit.parameters(dynamics);
                Symbolics.scalarize.(collect(Iterators.flatmap(c->Iterators.flatmap(Symbolics.get_variables, c.eqn), constraints)));
                Symbolics.get_variables(nonconvex_constraints);
                Symbolics.get_variables(terminal_cost);
                vcat((Symbolics.get_variables.(values(bcs)))...)])
        states = filter(Base.Fix2(isdynamics, false), vars)
        controls = filter(Base.Fix2(iscontrols, false), vars)
        parameters = filter(Base.Fix2(istunable, false), vars)
        @assert length(vars) == length(union(states, controls, parameters)) "Some variables have not been tagged!"
        
        tc = Symbolics.unwrap(terminal_cost)
        println(tc)
        println(terminal_cost)

        return new(states, 
            controls, 
            parameters, 
            unwrap_dict(scale_advice), 
            dynamics, 
            constraints, 
            tc, 
            running_cost, 
            Symbolics.unwrap.(nonconvex_constraints), 
            Dict(k => Symbolics.unwrap.(v) for (k,v) in bcs), 
            initalizer)
    end
end

function vartype(prob::SCPtProblem, p::SymbolicUtils.BasicSymbolic)
    if symin(p, prob.states) return :state 
    elseif symin(p, prob.controls) return :input
    elseif symin(p, prob.parameters) return :parameter
    else
        throw("Invalid variable for problem $p")
    end
end

#TODO: make the index uniquely identify which thingie it came from (e.g. StateIndex)
function varindex(prob::SCPtProblem, p::SymbolicUtils.BasicSymbolic)
    if symin(p, prob.states) return indexof(p, prob.states)
    elseif symin(p, prob.controls) return indexof(p, prob.controls)
    elseif symin(p, prob.parameters) return indexof(p, prob.parameters)
    else
        throw("Invalid variable for problem $p")
    end
end