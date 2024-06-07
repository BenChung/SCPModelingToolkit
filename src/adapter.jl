function mtk_add_conic_constraints!(pbm, prob::SCPtProblem, cstrs::Vector{ConicConstraint}, algo)
    state_constraints = ConicConstraint[]
    control_constraints = ConicConstraint[]
    combined_constraints = ConicConstraint[]
    for cs in cstrs
        if cs.eqn(0.0, 1) isa Vector
            vars = collect(Iterators.flatten(Symbolics.get_variables.(cs.eqn(0.0, 1))))
        else
            vars = Symbolics.get_variables(cs.eqn(0.0, 1))
        end
        if all(!isnothing, indexof.(vars, ( [prob.states; prob.parameters], )))
            push!(state_constraints, cs)
        elseif all(!isnothing, indexof.(vars, ([prob.controls; prob.parameters],)))
            push!(control_constraints, cs)
        else
            push!(combined_constraints, cs)
            #throw("constraint $(cs.name) with equation $(cs.eqn) has both state and control variables");
        end
    end
    # this is SUPER LAME but for totally unclear reasons SCPt doesn't support constraints that couple controls (including slacks!) and states
    # so we cheat! and keep a big array of states and constraint variables that are given to the X and U set setters

    states = Dict{Int, Any}()
    if length(state_constraints) > 0
        local symbolic_state_vars = [prob.states; prob.parameters]
        problem_set_X!(
            pbm,
            (t, k, x, p, pbm, ocp) -> begin 
                common = (pbm, ocp, algo)
                states[k] = x
                for cs in state_constraints
                    define_conic_constraint!(
                        common...,
                        cs.cone,
                        cs.name,
                        (x...,p...),
                        (args...) -> begin 
                            expr = Symbolics.value.(Base.Fix2(substitute, Dict(symbolic_state_vars .=> args)).(cs.eqn(t,k)))
                            if any(e-> e isa JuMP.AffExpr, expr)
                                expr = [e isa JuMP.AffExpr ? e : JuMP.AffExpr(e) for e in expr]
                            end
                            return expr
                        end,
                    )
                end
        end)
    if length(control_constraints) > 0
        local symbolic_control_vars = [prob.controls; prob.parameters]
        local combined_control_vars = [prob.states; prob.controls; prob.parameters]
        problem_set_U!(
            pbm,
            (t, k, u, p, pbm, ocp) -> begin 
                common = (pbm, ocp, algo)
                x = states[k]
                for cs in control_constraints
                    define_conic_constraint!(
                        common...,
                        cs.cone,
                        cs.name,
                        (u...,p...),
                        (args...) -> begin 
                            expr = Symbolics.value.(Base.Fix2(substitute, Dict(symbolic_control_vars .=> args)).(cs.eqn(t,k)))
                            if any(e-> e isa JuMP.AffExpr, expr)
                                expr = [e isa JuMP.AffExpr ? e : JuMP.AffExpr(e) for e in expr]
                            end
                            return expr
                        end,
                    )
                end
                for cs in combined_constraints
                    if value(x) isa Vector{Float64}
                        continue 
                    end
                    define_conic_constraint!(
                        common...,
                        cs.cone,
                        cs.name,
                        (x..., u..., p...),
                        (args...) -> Symbolics.value.(Base.Fix2(substitute, Dict(combined_control_vars .=> args)).(cs.eqn(t,k))),
                    )
                end
        end)
    end
    end
end

function mtk_add_nonconvex_constraints!(pbm, prob::SCPtProblem, sys, cstrs, algo)
    sfunc,_ = build_function(cstrs, prob.states, prob.controls, prob.parameters; expression=Val{false}())
    jfunc,_ = build_function(Symbolics.jacobian(cstrs, dynamics_variables(sys)), 
        prob.states, 
        prob.controls, 
        prob.parameters; expression=Val{false}())
    # Constraint s
    _q__s = (t, k, x, u, p, pbm) -> sfunc(x, u, p)
    # Jacobian ds/dx
    _q__C = (t, k, x, u, p, pbm) -> jfunc(x, u, p)
    
    if algo == :scvx
        problem_set_s!(pbm, algo, _q__s, _q__C)
    else
        throw("non-scvx not supported yet")
        #_q___s = (t, k, x, p, pbm) -> _q__s(t, k, x, nothing, p, pbm)
        #_q___C = (t, k, x, p, pbm) -> _q__C(t, k, x, nothing, p, pbm)
        #problem_set_s!(pbm, algo, _q___s, _q___C)
    end
end
function mtk_set_bc!(pbm, prob::SCPtProblem, kind, bc)
    g,_ = build_function(bc, prob.states, prob.parameters; expression=Val{false}())
    jfunc,_ = build_function(Symbolics.jacobian(bc, prob.states), 
                            prob.states, 
                            prob.parameters; expression=Val{false}())
    jpfunc,_ = build_function(Symbolics.jacobian(bc, prob.parameters), 
                            prob.states, 
                            prob.parameters; expression=Val{false}())


    problem_set_bc!(
        pbm,
        kind,
        # Constraint g
        (x, p, pbm) -> g(x, p),
        # Jacobian dg/dx
        (x, p, pbm) -> jfunc(x, p),
        # Jacobian dg/dp
        (x, p, pbm) -> jpfunc(x, p),
    )
end

function solve!(prob::SCPtProblem, solver; 
        parameters=Dict(),
        modify=nothing,
        N = 30, 
        Nsub = 15, 
        iter_max = 50,
        trials = 10,
        wvc = 1e3,
        wtr = 0.1,
        ϵ_abs = 1e-5,
        ϵ_rel = 0.01/100,
        feas_tol = 1e-3,
        q_tr = 8,
        q_exit = 2,
        solver_options = Dict("verbose" => false),
        use_forwarddiff=false)
    algo = :scvx
    pbm = TrajectoryProblem(prob)

    problem_set_dims!(pbm, length(prob.states), length(prob.controls), length(prob.parameters))

    if length(prob.scale_advice) > 0
        for (var, (min, max)) in prob.scale_advice
            problem_advise_scale!(pbm, vartype(prob, var), varindex(prob, var), (min, max))
        end
    end
    
    problem_set_terminal_cost!(pbm, build_function(prob.terminal_cost, prob.states, prob.parameters, expression=Val{false}))
    problem_set_running_cost!(pbm, algo, (t, k, x, u, p, pbm) -> Symbolics.value(substitute(prob.running_cost(t,k), [prob.states .=> x; collect(prob.controls .=> u); prob.parameters .=> p])))

    paramdict = merge(ModelingToolkit.defaults(prob.dynamics), parameters)
    ps = ModelingToolkit.varmap_to_vars(paramdict, 
        setdiff(setdiff(ModelingToolkit.parameters(prob.dynamics), ModelingToolkit.tunable_parameters(prob.dynamics)), prob.controls))

    func = first(internal_generate_function(prob.dynamics; expression=Val{false}))
    if !use_forwarddiff
        jac = first(mtk_jacobian(prob.dynamics, prob.states; expression=Val{false}))
        control_jac = first(mtk_jacobian(prob.dynamics, prob.controls; expression=Val{false}))
        param_jac = first(mtk_jacobian(prob.dynamics, prob.parameters; expression=Val{false}))
    else
        jac = (x,u,p,t) -> begin 
            sens = ForwardDiff.jacobian((x) -> func(x, u, p, ps, t), x)
            return sens
        end
        control_jac = (x,u,p,t) -> ForwardDiff.jacobian((u) -> func(x, u, p, ps, t), u)
        param_jac = (x,u,p,t) -> ForwardDiff.jacobian((p) -> func(x, u, p, ps, t), p) 
    end
    problem_set_dynamics!(
        pbm,
        # Dynamics f
        (t, k, x, u, p, pbm) -> begin
            res = func(x, u, p, ps, t)
            return res
        end,
        # Jacobian df/dx
        (t, k, x, u, p, pbm) -> begin
            res = jac(x, u, p, t)
            return res
        end,
        # Jacobian df/du
        (t, k, x, u, p, pbm) -> begin
            res = control_jac(x, u, p, t)
            return res
            end,
        # Jacobian df/dp
        (t, k, x, u, p, pbm) -> begin 
            res = param_jac(x, u, p, t)
            return res
        end,
    )

    mtk_add_conic_constraints!(pbm, prob, prob.constraints, algo)
    if !isempty(prob.nonconvex_constraints)
        mtk_add_nonconvex_constraints!(pbm, prob, prob.dynamics, prob.nonconvex_constraints, algo)
    end
    for (type, cst) in prob.bcs
        mtk_set_bc!(pbm, prob, type, cst)
    end
    if !isnothing(modify)
        modify(pbm)
    end
    initialize(pbm, prob, prob.initalizer)
    disc_method = FOH
    pars = PTR.Parameters(
        N,
        Nsub,
        iter_max,
        disc_method,
        wvc,
        wtr,
        ϵ_abs,
        ϵ_rel,
        feas_tol,
        q_tr,
        q_exit,
        solver,
        solver_options,
    )
    pbm.callback! = prob.callback
    pbm = PTR.create(pars, pbm);
    
    return PTR.solve(pbm)
end
