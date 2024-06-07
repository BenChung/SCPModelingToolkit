
struct ConicConstraint
    name::String
    eqn::Function # (t, k) -> Vector{SymbolicUtils.BasicSymbolic}
    cone::SCPToolbox.Parser.ConicLinearProgram.SupportedCone
    ConicConstraint(name::String, eqn::SymbolicUtils.BasicSymbolic, cone::SCPToolbox.Parser.ConicLinearProgram.SupportedCone) = new(name, (t, k) -> [eqn], cone)
    ConicConstraint(name::String, eqn::Num, cone::SCPToolbox.Parser.ConicLinearProgram.SupportedCone) = new(name, (t, k) -> [Symbolics.unwrap(eqn)], cone)
    ConicConstraint(name::String, eqns::Vector{SymbolicUtils.BasicSymbolic}, cone::SCPToolbox.Parser.ConicLinearProgram.SupportedCone) = new(name, (t, k) -> eqns, cone)
    ConicConstraint(name::String, eqns::Vector{Num}, cone::SCPToolbox.Parser.ConicLinearProgram.SupportedCone) = new(name, (t, k) -> Symbolics.unwrap.(eqns), cone)
    ConicConstraint(name::String, eqn_generator::Function, cone::SCPToolbox.Parser.ConicLinearProgram.SupportedCone) = new(name, eqn_generator, cone)
end
abstract type Initialization end

struct SCPtProblem{F}
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
    callback::F
    function SCPtProblem(; 
        scale_advice::Dict{Num, Tuple{Float64, Float64}} = Dict(), 
        dynamics::ModelingToolkit.ODESystem,
        constraints::Vector{ConicConstraint},
        terminal_cost::Union{Float64, Num},
        running_cost::Function, 
        nonconvex_constraints::Vector{Num},
        bcs::Dict{Symbol, Vector{Num}},
        initalizer::Initialization,
        slacks::Vector = [],
        callback::F = x->nothing) where F
        vars = unique([ModelingToolkit.states(dynamics);
                SCPModelingToolkit.control_variables(dynamics);
                filter(Base.Fix2(istunable, false), ModelingToolkit.parameters(dynamics));
                Symbolics.scalarize.(collect(Iterators.flatmap(c->Iterators.flatmap(Symbolics.get_variables, c.eqn(0.0, 1)), constraints)));
                Symbolics.get_variables(nonconvex_constraints);
                Symbolics.get_variables(terminal_cost);
                vcat((Symbolics.get_variables.(values(bcs)))...);
                reduce(vcat, Symbolics.get_variables.(slacks); init=[])])
        states = filter(Base.Fix2(isdynamics, false), vars)
        controls = filter(Base.Fix2(iscontrols, false), vars)
        parameters = filter(Base.Fix2(istunable, false), vars)
        if length(vars) != length(union(states, controls, parameters)) 
            @debug "Some variables have not been tagged! \n states: $states \n controls: $controls \n parameters: $parameters \n vars: $vars"
        end
        
        tc = Symbolics.unwrap(terminal_cost)
        return new{F}(states, 
            controls, 
            parameters, 
            unwrap_dict(scale_advice), 
            dynamics, 
            constraints, 
            tc, 
            running_cost, 
            Symbolics.unwrap.(nonconvex_constraints), 
            Dict(k => Symbolics.unwrap.(v) for (k,v) in bcs), 
            initalizer,
            callback)
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