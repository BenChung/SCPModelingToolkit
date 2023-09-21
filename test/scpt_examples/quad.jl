u_max = 23.2
u_min = 0.6
tilt_max = deg2rad(60)

# >> Environment <<
g = [0, 0, -9.81]
obs = [
    Ellipsoid(diagm([2.0; 2.0; 0.0]), [1.0; 2.0; 0.0]),
    Ellipsoid(diagm([1.5; 1.5; 0.0]), [2.0; 5.0; 0.0]),
]

# >> Trajectory <<
r0 = zeros(3)
rf = zeros(3)
rf[1:2] = [2.5; 6.0]
v0 = zeros(3)
vf = zeros(3)
tf_min = 0.0
tf_max = 2.5
γ = 0.0


ModelingToolkit.@variables t x(t)[1:3] [dynamics=true] v(t)[1:3] [dynamics=true]
ModelingToolkit.@variables u(t)[1:3] [scpcontrol = true] σ(t) [scpcontrol = true]
@parameters τ [tunable = true]
D = Differential(t)

acceleration = D.(v) .~ ([0,0,-9.8] + u)*τ
velocity = D.(x) .~ v*τ

@named sys = ODESystem(Symbolics.scalarize.([velocity; acceleration]), t, [x; v; u], [τ])
sys,_ = structural_simplify(sys, (SCPModelingToolkit.control_variables(sys), SCPModelingToolkit.dynamics_variables(sys)))

function (E::Ellipsoid)(r)
    y = norm(E.H * (r - E.c))
    return y
end

prob = SCPModelingToolkit.SCPtProblem(
    scale_advice = Dict(τ => (tf_min, tf_min + 1.0 * (tf_max - tf_min))),
    dynamics = sys,
    constraints = SCPModelingToolkit.ConicConstraint[
        SCPModelingToolkit.ConicConstraint("min_accel", u_min - σ, NONPOS),
        SCPModelingToolkit.ConicConstraint("max_accel", σ - u_max, NONPOS),
        SCPModelingToolkit.ConicConstraint("lcvx_equality", [σ; u], SOC),
        SCPModelingToolkit.ConicConstraint("max_tilt", σ * cos(tilt_max) - u[3], NONPOS),
        SCPModelingToolkit.ConicConstraint("max_duration", τ - tf_max, NONPOS),
        SCPModelingToolkit.ConicConstraint("min_duration", tf_min - τ, NONPOS),
    ],
    terminal_cost = γ * (τ/tf_max)^2,
    running_cost = (t,k) -> (1 - γ)*(σ/norm(g))^2,
    nonconvex_constraints = [1 - obs(Symbolics.scalarize(x)) for obs in obs],
    bcs = Dict{Symbol, Vector{Num}}(
        :ic => [x;v] - [r0; v0],
        :tc => [x;v] - [rf; vf]
    ),
    initalizer = SCPModelingToolkit.StraightlineInterpolate(
        Dict{Num, Float64}(Symbolics.scalarize.([x .=> r0; v .=> v0; u .=> -g; σ => norm(g)])), 
        Dict{Num, Float64}(Symbolics.scalarize.([x .=> rf; v .=> vf; u .=> -g; σ => norm(g)])),
        Dict{Num, Float64}([τ => 0.5 * (tf_min + tf_max)]))
)
SCPModelingToolkit.solve!(prob, ECOS)