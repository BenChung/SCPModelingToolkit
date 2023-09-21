
indexof(sym, syms) = findfirst(isequal(sym), syms)
symin(sym, syms) = any(isequal(sym), syms)
function mtk_jacobian(sys::ModelingToolkit.AbstractODESystem, dvs = states(sys), ps = ModelingToolkit.tunable_parameters(sys);
    simplify = false, sparse = false, kwargs...)
    jac = calculate_jacobian(sys;dvs=dvs, simplify = simplify, sparse = sparse)
    pre = ModelingToolkit.get_preprocess_constants(jac)
    return build_function(jac, ModelingToolkit.outputs(sys), ModelingToolkit.inputs(sys), ps, ModelingToolkit.get_iv(sys); postprocess_fbody = pre, kwargs...)
end


function internal_generate_function(sys::ModelingToolkit.ODESystem, ps = ModelingToolkit.tunable_parameters(sys);
    has_difference = false,
    kwargs...)
    eqs = [eq for eq in equations(sys) if !ModelingToolkit.isdifferenceeq(eq)]
    ModelingToolkit.check_operator_variables(eqs, Differential)
    ModelingToolkit.check_lhs(eqs, Differential, Set(ModelingToolkit.states(sys)))
    # substitute x(t) by just x
    reorder_eqns = indexof.([arguments(eq.lhs)[1] for eq in sys.eqs if ModelingToolkit.isdifferential(eq.lhs)], (states(sys), ))
    rhss = [eq.rhs for eq in eqs][reorder_eqns]

    # TODO: add an optional check on the ordering of observed equations
    x = map(x -> ModelingToolkit.time_varying_as_func(ModelingToolkit.value(x), sys), ModelingToolkit.outputs(sys))
    u = map(x -> ModelingToolkit.time_varying_as_func(ModelingToolkit.value(x), sys), ModelingToolkit.inputs(sys))
    p = map(x -> ModelingToolkit.time_varying_as_func(ModelingToolkit.value(x), sys), ModelingToolkit.tunable_parameters(sys))
    t = ModelingToolkit.get_iv(sys)

    pre, sol_states = ModelingToolkit.get_substitutions_and_solved_states(sys, no_postprocess = has_difference)

    build_function(rhss, x, u, p, t; postprocess_fbody = pre, states = sol_states,
        kwargs...)
end

unwrap_dict(d) = Dict(Symbolics.unwrap(k) => v for (k,v) in d)