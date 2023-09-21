
struct DynamicVariables end
Symbolics.option_to_metadata_type(::Val{:dynamics}) = DynamicVariables
struct SCPControlVariable end
Symbolics.option_to_metadata_type(::Val{:scpcontrol}) = SCPControlVariable

function isdynamics(x, default = false)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, DynamicVariables, default)
end
isdynamics(x::Num, args...) = isdynamics(Symbolics.unwrap(x), args...)
function dynamics_variables(sys, vs = states(sys); default = false)
    filter(x -> isdynamics(x, default), vs)
end

function iscontrols(x, default = false)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, SCPControlVariable, default)
end
iscontrols(x::Num, args...) = iscontrols(Symbolics.unwrap(x), args...)
function control_variables(sys, vs = [states(sys); ModelingToolkit.parameters(sys)]; default = false)
    filter(x -> iscontrols(x, default), vs)
end