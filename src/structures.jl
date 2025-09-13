################ Introduce the model structures ###############

"""
Parent model type.
"""
abstract type DSGEModel end

abstract type REModel <: DSGEModel end
abstract type REModelPartial <: DSGEModel end

"""
Parent type for the primative or key components of a DSGEModel.
"""
abstract type DSGEModelPrimatives end

"""
Struct containing the primatives for a rational expectations model.
"""
struct REModelPrimatives{Q<:AbstractString} <: DSGEModelPrimatives

    # This structure contains the crucial model information extracted from a
    # model file.

    states::Array{Q,1}
    jumps::Array{Q,1}
    shocks::Array{Q,1}
    variables::Array{Q,1}
    parameters::Array{Q,1}
    parametervalues::Array{Q,1}
    equations::Array{Q,1}
    unassigned_parameters::Array{Q,1}
    solvers::Q

end

"""
Struct for a fully specified linear rational expectations model, a subtype of REModel.
    
Contains the model information in a form that the solvers can work with.
"""
struct REModelLinear{S<:Integer,Q<:AbstractString} <: REModel

    number_states::S
    number_jumps::S
    number_shocks::S
    number_variables::S
    number_equations::S
    variables::OrderedDict{Q,S}
    nlsolve_static_function::Function
    static_function::Function
    dynamic_function::Function

end

"""
Struct for a fully specified rational expectations model, a subtype of REModel.
    
Contains the model information in a form that the solvers can work with.
"""
struct REModelPert{S<:Integer,Q<:AbstractString} <: REModel

    number_states::S
    number_jumps::S
    number_shocks::S
    number_variables::S
    number_equations::S
    variables::OrderedDict{Q,S}
    nlsolve_static_function::Function
    static_function::Function
    dynamic_function::Function
    each_eqn_function::Array{Function,1}

end

"""
Struct for a fully specified rational expectations model, a subtype of REModel.
    
Contains the model information in a form that the solvers can work with.
"""
struct REModelProj{S<:Integer,Q<:AbstractString} <: REModel

    number_states::S
    number_jumps::S
    number_shocks::S
    number_variables::S
    number_equations::S
    jumps_approximated::Array{S,1}
    eqns_approximated::Array{S,1}
    derivs_approximated_num::Array{S,1}
    derivs_approximated_den::Array{S,1}
    eqns_with_derivs::Array{S,1}
    variables::OrderedDict{Q,S}
    nlsolve_static_function::Function
    static_function::Function
    dynamic_function::Function
    each_eqn_function::Array{Function,1}
    closure_function_chebyshev::Function
    closure_function_smolyak::Function
    closure_function_hcross::Function
    closure_function_piecewise::Function

end

"""
Struct for a fully specified rational expectations model, a subtype of REModel.
    
Contains the model information in a form that the solvers can work with.
"""
struct REModelAny{S<:Integer,Q<:AbstractString} <: REModel

    number_states::S
    number_jumps::S
    number_shocks::S
    number_variables::S
    number_equations::S
    jumps_approximated::Array{S,1}
    eqns_approximated::Array{S,1}
    derivs_approximated_num::Array{S,1}
    derivs_approximated_den::Array{S,1}
    eqns_with_derivs::Array{S,1}
    variables::OrderedDict{Q,S}
    nlsolve_static_function::Function
    static_function::Function
    dynamic_function::Function
    each_eqn_function::Array{Function,1}
    closure_function_chebyshev::Function
    closure_function_smolyak::Function
    closure_function_hcross::Function
    closure_function_piecewise::Function

end

"""
Struct for a partially specified linear rational expectations model, a subtype of REModel.

The model is partially specified in that one or more parameters lack values.  These values must be
assigned using the assign_parameters() function in order for the model to become fully specified.
"""
struct REModelPartialLinear{S<:Integer,Q<:AbstractString} <: REModelPartial

    number_states::S
    number_jumps::S
    number_shocks::S
    number_variables::S
    number_equations::S
    variables::OrderedDict{Q,S}
    nlsolve_static_function::Function
    static_function::Function
    dynamic_function::Function
    unassigned_parameters::Array{Q,1}

end

"""
Struct for a partially specified rational expectations model, a subtype of REModel.

The model is partially specified in that one or more parameters lack values.  These values must be
assigned using the assign_parameters() function in order for the model to become fully specified.
"""
struct REModelPartialPert{S<:Integer,Q<:AbstractString} <: REModelPartial

    number_states::S
    number_jumps::S
    number_shocks::S
    number_variables::S
    number_equations::S
    variables::OrderedDict{Q,S}
    nlsolve_static_function::Function
    static_function::Function
    dynamic_function::Function
    each_eqn_function::Array{Function,1}
    unassigned_parameters::Array{Q,1}
    
end

"""
Struct for a partially specified rational expectations model, a subtype of REModel.

The model is partially specified in that one or more parameters lack values.  These values must be
assigned using the assign_parameters() function in order for the model to become fully specified.
"""
struct REModelPartialProj{S<:Integer,Q<:AbstractString} <: REModelPartial

    number_states::S
    number_jumps::S
    number_shocks::S
    number_variables::S
    number_equations::S
    jumps_approximated::Array{S,1}
    eqns_approximated::Array{S,1}
    derivs_approximated_num::Array{S,1}
    derivs_approximated_den::Array{S,1}
    eqns_with_derivs::Array{S,1}
    variables::OrderedDict{Q,S}
    nlsolve_static_function::Function
    static_function::Function
    dynamic_function::Function
    each_eqn_function::Array{Function,1}
    closure_function_chebyshev::Function
    closure_function_smolyak::Function
    closure_function_hcross::Function
    closure_function_piecewise::Function
    unassigned_parameters::Array{Q,1}

end

"""
Struct for a partially specified rational expectations model, a subtype of REModel.

The model is partially specified in that one or more parameters lack values.  These values must be
assigned using the assign_parameters() function in order for the model to become fully specified.
"""
struct REModelPartialAny{S<:Integer,Q<:AbstractString} <: REModelPartial

    number_states::S
    number_jumps::S
    number_shocks::S
    number_variables::S
    number_equations::S
    jumps_approximated::Array{S,1}
    eqns_approximated::Array{S,1}
    derivs_approximated_num::Array{S,1}
    derivs_approximated_den::Array{S,1}
    eqns_with_derivs::Array{S,1}
    variables::OrderedDict{Q,S}
    nlsolve_static_function::Function
    static_function::Function
    dynamic_function::Function
    each_eqn_function::Array{Function,1}
    closure_function_chebyshev::Function
    closure_function_smolyak::Function
    closure_function_hcross::Function
    closure_function_piecewise::Function
    unassigned_parameters::Array{Q,1}

end

struct REModelFullDetail{S<:Integer,Q<:AbstractString} <: DSGEModelPrimatives

    number_states::S
    number_jumps::S
    number_shocks::S
    number_variables::S
    number_equations::S
    jumps_approximated::Array{S,1}
    eqns_approximated::Array{S,1}
    derivs_approximated_num::Array{S,1}
    derivs_approximated_den::Array{S,1}
    eqns_with_derivs::Array{S,1}
    variables::OrderedDict{Q,S}
    nlsolve_static_function::Function
    static_function::Function
    dynamic_function::Function
    each_eqn_function::Array{Function,1}
    closure_function_chebyshev::Function
    closure_function_smolyak::Function
    closure_function_hcross::Function
    closure_function_piecewise::Function
    unassigned_parameters::Array{Q,1}
    solvers::Q

end

"""
Makes a copy of any DSGEModel sub-type (REModel and REModelPartial).

```
new_model = copy_model_struct(model)
```
"""
function copy_model_struct(mod::REModelLinear)

    new_model = REModelLinear(mod.number_states,mod.number_jumps,mod.number_shocks,mod.number_variables,mod.number_equations,mod.variables,mod.nlsolve_static_function,
                        mod.static_function,mod.dynamic_function)

    return new_model

end

function copy_model_struct(mod::REModelPartialLinear)

    new_model = REModelPartialLinear(mod.number_states,mod.number_jumps,mod.number_shocks,mod.number_variables,mod.number_equations,mod.variables,mod.nlsolve_static_function,
                               mod.static_function,mod.dynamic_function,mod.unassigned_parameters)

    return new_model

end

function copy_model_struct(mod::REModelPert)

    new_model = REModelPert(mod.number_states,mod.number_jumps,mod.number_shocks,mod.number_variables,mod.number_equations,mod.variables,mod.nlsolve_static_function,
                        mod.static_function,mod.dynamic_function,mod.each_eqn_function)

    return new_model

end

function copy_model_struct(mod::REModelPartialPert)

    new_model = REModelPartialPert(mod.number_states,mod.number_jumps,mod.number_shocks,mod.number_variables,mod.number_equations,mod.variables,mod.nlsolve_static_function,
                               mod.static_function,mod.dynamic_function,mod.each_eqn_function,mod.unassigned_parameters)

    return new_model

end

function copy_model_struct(mod::REModelProj)

    new_model = REModelProj(mod.number_states,mod.number_jumps,mod.number_shocks,mod.number_variables,mod.number_equations,mod.jumps_approximated,mod.eqns_approximated,mod.derivs_approximated_num,mod.derivs_approximated_den,mod.eqns_with_derivs,mod.variables,mod.nlsolve_static_function,
                        mod.static_function,mod.dynamic_function,mod.each_eqn_function,mod.closure_function_chebyshev,mod.closure_function_smolyak,mod.closure_function_hcross,mod.closure_function_piecewise)

    return new_model

end

function copy_model_struct(mod::REModelPartialProj)

    new_model = REModelPartialProj(mod.number_states,mod.number_jumps,mod.number_shocks,mod.number_variables,mod.number_equations,mod.jumps_approximated,mod.eqns_approximated,mod.derivs_approximated_num,mod.derivs_approximated_den,mod.eqns_with_derivs,mod.variables,mod.nlsolve_static_function,
                               mod.static_function,mod.dynamic_function,mod.each_eqn_function,mod.closure_function_chebyshev,mod.closure_function_smolyak,mod.closure_function_hcross,mod.closure_function_piecewise,mod.unassigned_parameters)

    return new_model

end

function copy_model_struct(mod::REModelAny)

    new_model = REModelAny(mod.number_states,mod.number_jumps,mod.number_shocks,mod.number_variables,mod.number_equations,mod.jumps_approximated,mod.eqns_approximated,mod.derivs_approximated_num,mod.derivs_approximated_den,mod.eqns_with_derivs,mod.variables,mod.nlsolve_static_function,
                        mod.static_function,mod.dynamic_function,mod.each_eqn_function,mod.closure_function_chebyshev,mod.closure_function_smolyak,mod.closure_function_hcross,mod.closure_function_piecewise)

    return new_model

end

function copy_model_struct(mod::REModelPartialAny)

    new_model = REModelPartialAny(mod.number_states,mod.number_jumps,mod.number_shocks,mod.number_variables,mod.number_equations,mod.jumps_approximated,mod.eqns_approximated,mod.derivs_approximated_num,mod.derivs_approximated_den,mod.eqns_with_derivs,mod.variables,mod.nlsolve_static_function,
                               mod.static_function,mod.dynamic_function,mod.each_eqn_function,mod.closure_function_chebyshev,mod.closure_function_smolyak,mod.closure_function_hcross,mod.closure_function_piecewise,
                               mod.unassigned_parameters)

    return new_model

end

############### Introduce the solutionscheme structures #######################

abstract type SolutionScheme end
abstract type ProjectionScheme <: SolutionScheme end
abstract type ProjectionSchemeDet <: ProjectionScheme end
abstract type ProjectionSchemeStoch <: ProjectionScheme end
abstract type ProjectionSchemeOBCDet <: ProjectionScheme end
abstract type ProjectionSchemeOBCStoch <: ProjectionScheme end

struct PerturbationScheme{T<:Real,Q<:AbstractString} <: SolutionScheme

    steady_state::Array{T,1}
    cutoff::T
    order::Q

end

struct ChebyshevSchemeDet{T<:Real,S<:Integer} <: ProjectionSchemeDet

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    node_number::Union{S,Array{S,1}}
    order::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct ChebyshevSchemeOBCDet{T<:Real,S<:Integer} <: ProjectionSchemeOBCDet

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    node_number::Union{S,Array{S,1}}
    order::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    lb::Array{T,1}
    ub::Array{T,1}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct ChebyshevSchemeStoch{T<:Real,S<:Integer} <: ProjectionSchemeStoch

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    node_number::Union{S,Array{S,1}}
    num_quad_nodes::S
    order::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct ChebyshevSchemeOBCStoch{T<:Real,S<:Integer} <: ProjectionSchemeOBCStoch

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    node_number::Union{S,Array{S,1}}
    num_quad_nodes::S
    order::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    lb::Array{T,1}
    ub::Array{T,1}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct SmolyakSchemeDet{T<:Real,S<:Integer} <: ProjectionSchemeDet

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    layer::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct SmolyakSchemeOBCDet{T<:Real,S<:Integer} <: ProjectionSchemeOBCDet

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    layer::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    lb::Array{T,1}
    ub::Array{T,1}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct SmolyakSchemeStoch{T<:Real,S<:Integer} <: ProjectionSchemeStoch

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    num_quad_nodes::S
    layer::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct SmolyakSchemeOBCStoch{T<:Real,S<:Integer} <: ProjectionSchemeOBCStoch

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    num_quad_nodes::S
    layer::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    lb::Array{T,1}
    ub::Array{T,1}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct HyperbolicCrossSchemeDet{T<:Real,S<:Integer} <: ProjectionSchemeDet

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    layer::S
    n::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct HyperbolicCrossSchemeOBCDet{T<:Real,S<:Integer} <: ProjectionSchemeOBCDet

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    layer::S
    n::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    lb::Array{T,1}
    ub::Array{T,1}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct HyperbolicCrossSchemeStoch{T<:Real,S<:Integer} <: ProjectionSchemeStoch

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    num_quad_nodes::S
    layer::S
    n::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct HyperbolicCrossSchemeOBCStoch{T<:Real,S<:Integer} <: ProjectionSchemeOBCStoch

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    num_quad_nodes::S
    layer::S
    n::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    lb::Array{T,1}
    ub::Array{T,1}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct PiecewiseLinearSchemeDet{T<:Real,S<:Integer} <: ProjectionSchemeDet

    initial_guess::Union{T,Array{T,1}}
    node_number::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct PiecewiseLinearSchemeOBCDet{T<:Real,S<:Integer} <: ProjectionSchemeOBCDet

    initial_guess::Union{T,Array{T,1}}
    node_number::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    lb::Array{T,1}
    ub::Array{T,1}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct PiecewiseLinearSchemeStoch{T<:Real,S<:Integer} <: ProjectionSchemeStoch

    initial_guess::Union{T,Array{T,1}}
    node_number::Union{S,Array{S,1}}
    num_quad_nodes::S
    domain::Union{Array{T,1},Array{T,2}}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

struct PiecewiseLinearSchemeOBCStoch{T<:Real,S<:Integer} <: ProjectionSchemeOBCStoch

    initial_guess::Union{T,Array{T,1}}
    node_number::Union{S,Array{S,1}}
    num_quad_nodes::S
    domain::Union{Array{T,1},Array{T,2}}
    lb::Array{T,1}
    ub::Array{T,1}
    ftol::T
    xtol::T
    maxiters::S
    method::Symbol

end

##################### Introduce the solution structures ########################

abstract type ModelSolution end

abstract type PerturbationSolution <: ModelSolution end
abstract type PerturbationSolutionDet <: PerturbationSolution end
abstract type PerturbationSolutionStoch <: PerturbationSolution end

abstract type ProjectionSolution <: ModelSolution end
abstract type ProjectionSolutionDet <: ProjectionSolution end
abstract type ProjectionSolutionStoch <: ProjectionSolution end

struct FirstOrderSolutionDet{T<:Real,S<:Integer} <: PerturbationSolutionDet

    # x(t+1) = hx*x(t)
    #   y(t) = gx*x(t)

    hbar::Union{T,Array{T,1}}           # steady state values for predetermined variables
    hx::Union{Array{T,2},Array{T,1}}    # Transition matrix for predetermined variables
    gbar::Union{T,Array{T,1}}           # steady state values for nonpredetermined variables
    gx::Union{Array{T,2},Array{T,1}}    # Decision rule matrix linking nonpredetermined variables to predetermined variables
    grc::S                              # Number of eigenvalues greater than cutoff
    soln_type::String                   # "determinate", "indeterminate", or "explosive"

end

struct FirstOrderSolutionStoch{T<:Real,S<:Integer} <: PerturbationSolutionStoch

    # x(t+1) = hx*x(t) + k*v(t+1)
    #   y(t) = gx*x(t)

    hbar::Union{T,Array{T,1}}           # steady state values for predetermined variables
    hx::Union{Array{T,2},Array{T,1}}    # Transition matrix for predetermined variables
    k::Union{Array{T,2},Array{T,1}}     # Innovation loading matrix
    gbar::Union{T,Array{T,1}}           # steady state values for nonpredetermined variables
    gx::Union{Array{T,2},Array{T,1}}    # Decision rule matrix linking nonpredetermined variables to predetermined variables
    sigma::Union{Array{T,2},Array{T,1}} # Innovation variance-covariance matrix
    grc::S                              # Number of eigenvalues greater than cutoff
    soln_type::String                   # "determinate", "indeterminate", or "explosive"

end

struct SecondOrderSolutionDet{T<:Real,S<:Integer} <: PerturbationSolutionDet

    # x(t+1) = hx*x(t) + (1/2)*[kron(I,x(t))]'hxx*[kron(I,x(t))]
    #   y(t) = gx*x(t) + (1/2)*[kron(I,x(t))]'gxx*[kron(I,x(t))]

    hbar::Union{T,Array{T,1}}           # steady state values for predetermined variables
    hx::Union{Array{T,2},Array{T,1}}    # Linear component in predetermined block
    hxx::Array{T,2}                     # Quadratic component in predetermined block
    gbar::Union{T,Array{T,1}}           # steady state values for nonpredetermined variables
    gx::Union{Array{T,2},Array{T,1}}    # Linear component in non-predetermined block
    gxx::Array{T,2}                     # Quadratic component in non-predetermined block
    grc::S                              # Number of eigenvalues greater than cutoff
    soln_type::String                   # "determinate", "indeterminate", or "explosive"

end

struct SecondOrderSolutionStoch{T<:Real,S<:Integer} <: PerturbationSolutionStoch

    # x(t+1) = hx*x(t) + (1/2)*hss + (1/2)*[kron(I,x(t))]'hxx*[kron(I,x(t))] + k*v(t+1)
    #   y(t) = gx*x(t) + (1/2)*gss + (1/2)*[kron(I,x(t))]'gxx*[kron(I,x(t))]

    hbar::Union{T,Array{T,1}}           # steady state values for predetermined variables
    hx::Union{Array{T,2},Array{T,1}}    # Linear component in predetermined block
    hss::Array{T,1}                     # Intercepts in predetermined block
    hxx::Array{T,2}                     # Quadratic component in predetermined block
    k::Union{Array{T,2},Array{T,1}}     # Innovation loading matrix
    gbar::Union{T,Array{T,1}}           # steady state values for nonpredetermined variables
    gx::Union{Array{T,2},Array{T,1}}    # Linear component in non-predetermined block
    gss::Array{T,1}                     # Intercepts in predetermined block
    gxx::Array{T,2}                     # Quadratic component in non-predetermined block
    sigma::Union{Array{T,2},Array{T,1}} # Innovation variance-covariance matrix
    grc::S                              # Number of eigenvalues greater than cutoff
    soln_type::String                   # "determinate", "indeterminate", or "explosive"

end

struct ThirdOrderSolutionDet{T<:Real,S<:Integer} <: PerturbationSolutionDet

    # x(t+1) = hx*x(t) + (1/2)*hxx*[kron(x(t),x(t)] + (1/6)*hxxx*[kron(x(t),x(t),x(t))]

    #   y(t) = gx*x(t) + (1/2)*gxx*[kron(x(t),x(t)] + (1/6)*gxxx*[kron(x(t),x(t),x(t))]

    hbar::Union{T,Array{T,1}}              # steady state values for predetermined variables
    hx::Union{Array{T,2},Array{T,1}}       # Linear component in predetermined block
    hxx::Array{T,2}                        # Quadratic component in predetermined block
    hxxx::Array{T,2}                       # Third-order component in predetermined block
    gbar::Union{T,Array{T,1}}              # steady state values for nonpredetermined variables
    gx::Union{Array{T,2},Array{T,1}}       # Linear component in non-predetermined block
    gxx::Array{T,2}                        # Quadratic component in non-predetermined block
    gxxx::Array{T,2}                       # Third-order component in non-predetermined block
    grc::S                                 # Number of eigenvalues greater than cutoff
    soln_type::String                      # "determinate", "indeterminate", or "explosive"

end

struct ThirdOrderSolutionStoch{T<:Real,S<:Integer} <: PerturbationSolutionStoch

    # x(t+1) = hx*x(t) + (1/2)*hss + (1/2)*hxx*[kron(x(t),x(t)]
    #        + (1/6)*hsss + (3/6)*hssx*[x(t)] + (1/6)*hxxx*[kron(x(t),x(t),x(t))]
    #        + k*v(t+1)

    #   y(t) = gx*x(t) + (1/2)*gss + (1/2)*gxx*[kron(x(t),x(t)]
    #        +(1/6)*gsss + (3/6)*gssx*[x(t)] + (1/6)*gxxx*[kron(x(t),x(t),x(t))]

    hbar::Union{T,Array{T,1}}              # steady state values for predetermined variables
    hx::Union{Array{T,2},Array{T,1}}       # Linear component in predetermined block
    hss::Array{T,1}                        # Intercepts in predetermined block
    hxx::Array{T,2}                        # Quadratic component in predetermined block
    hsss::Array{T,1}                       # Third-order component in predetermined block
    hssx::Array{T,2}                       # Third-order component in predetermined block
    hxxx::Array{T,2}                       # Third-order component in predetermined block
    k::Union{Array{T,2},Array{T,1}}        # Innovation loading matrix
    gbar::Union{T,Array{T,1}}              # steady state values for nonpredetermined variables
    gx::Union{Array{T,2},Array{T,1}}       # Linear component in non-predetermined block
    gss::Array{T,1}                        # Intercepts in predetermined block
    gxx::Array{T,2}                        # Quadratic component in non-predetermined block
    gsss::Array{T,1}                       # Third-order component in non-predetermined block
    gssx::Array{T,2}                       # Third-order component in non-predetermined block
    gxxx::Array{T,2}                       # Third-order component in non-predetermined block
    sigma::Union{Array{T,2},Array{T,1}}    # Innovation variance-covariance matrix
    skewness::Union{Array{T,1},Array{T,2}} # Innovation skewness matrix
    grc::S                                 # Number of eigenvalues greater than cutoff
    soln_type::String                      # "determinate", "indeterminate", or "explosive"

end

struct FourthOrderSolutionDet{T<:Real,S<:Integer} <: PerturbationSolutionDet

    # x(t+1) = hx*x(t) + (1/2)*hxx*[kron(x(t),x(t)] + (1/6)*hxxx*[kron(x(t),x(t),x(t))] + (1/24)*hxxxx*[kron(x(t),x(t),x(t),x(t))]

    #   y(t) = gx*x(t) + (1/2)*gxx*[kron(x(t),x(t)] + (1/6)*gxxx*[kron(x(t),x(t),x(t))] + (1/24)*gxxxx*[kron(x(t),x(t),x(t),x(t))]

    hbar::Union{T,Array{T,1}}              # steady state values for predetermined variables
    hx::Union{Array{T,2},Array{T,1}}       # Linear component in predetermined block
    hxx::Array{T,2}                        # Quadratic component in predetermined block
    hxxx::Array{T,2}                       # Third-order component in predetermined block
    hxxxx::Array{T,2}                      # Fourth-order component in predetermined block
    gbar::Union{T,Array{T,1}}              # steady state values for nonpredetermined variables
    gx::Union{Array{T,2},Array{T,1}}       # Linear component in non-predetermined block
    gxx::Array{T,2}                        # Quadratic component in non-predetermined block
    gxxx::Array{T,2}                       # Third-order component in non-predetermined block
    gxxxx::Array{T,2}                      # Fourth-order component in non-predetermined block
    grc::S                                 # Number of eigenvalues greater than cutoff
    soln_type::String                      # "determinate", "indeterminate", or "explosive"

end

struct FourthOrderSolutionStoch{T<:Real,S<:Integer} <: PerturbationSolutionStoch

    # x(t+1) = hx*x(t) + (1/2)*hss + (1/2)*hxx*[kron(x(t),x(t)]
    #        + (3/6)*hssx*[x(t)] + (1/6)*hxxx*[kron(x(t),x(t),x(t))] 
    #        + (1/24)*hssss + (6/24)*hssxx*[kron(x(t),x(t))] + (1/24)*hxxxx*[kron(x(t),x(t),x(t),x(t))]
    #        + k*v(t+1)

    #   y(t) = gx*x(t) + (1/2)*gss + (1/2)*gxx*[kron(x(t),x(t)]
    #        + (3/6)*gssx*[x(t)] + (1/6)*gxxx*[kron(x(t),x(t),x(t))]
    #        + (1/24)*gssss + (6/24)*gssxx*[kron(x(t),x(t))] + (1/24)*gxxxx*[kron(x(t),x(t),x(t),x(t))]

    hbar::Union{T,Array{T,1}}              # steady state values for predetermined variables
    hx::Union{Array{T,2},Array{T,1}}       # Linear component in predetermined block
    hss::Array{T,1}                        # Intercepts in predetermined block
    hxx::Array{T,2}                        # Quadratic component in predetermined block
    hssx::Array{T,2}                       # Third-order component in predetermined block
    hxxx::Array{T,2}                       # Third-order component in predetermined block
    hssss::Array{T,1}                      # Fourth-order component in predetermined block
    hssxx::Array{T,2}                      # Fourth-order component in predetermined block
    hxxxx::Array{T,2}                      # Fourth-order component in predetermined block
    k::Union{Array{T,2},Array{T,1}}        # Innovation loading matrix
    gbar::Union{T,Array{T,1}}              # steady state values for nonpredetermined variables
    gx::Union{Array{T,2},Array{T,1}}       # Linear component in non-predetermined block
    gss::Array{T,1}                        # Intercepts in predetermined block
    gxx::Array{T,2}                        # Quadratic component in non-predetermined block
    gssx::Array{T,2}                       # Third-order component in non-predetermined block
    gxxx::Array{T,2}                       # Third-order component in non-predetermined block
    gssss::Array{T,1}                      # Fourth-order component in non-predetermined block
    gssxx::Array{T,2}                      # Fourth-order component in non-predetermined block
    gxxxx::Array{T,2}                      # Fourth-order component in non-predetermined block
    sigma::Union{Array{T,2},Array{T,1}}    # Innovation variance-covariance matrix
    grc::S                                 # Number of eigenvalues greater than cutoff
    soln_type::String                      # "determinate", "indeterminate", or "explosive"

end

struct ChebyshevSolutionDet{T<:Real,S<:Integer,N} <: ProjectionSolutionDet

    variables::Array{Array{T,N},1}       # Variables
    weights::Array{Array{T,N},1}         # Chebyshev polynomials
    nodes::Array{Array{T,1},1}           # Chebyshev nodes
    order::Union{S,Array{S,1}}           # Complete polynomial / tensor-product
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    iteration_count::S                   # Number of iterations needed for convergence
    node_generator::Function             # Function to generate the nodes

end

struct ChebyshevSolutionStoch{T<:Real,S<:Integer,N} <: ProjectionSolutionStoch

    variables::Array{Array{T,N},1}                   # Variables
    weights::Array{Array{T,N},1}                     # Chebyshev weights
    integrals::Union{Array{T},Array{Array{T,1},1}}   # Integrals for computing scaled weights
    nodes::Array{Array{T,1},1}                       # Chebyshev nodes
    order::Union{S,Array{S,1}}                       # Complete polynomial / tensor-product
    domain::Union{Array{T,2},Array{T,1}}             # Domain for state variables / state variable
    k::Union{Array{T,2},Array{T,1}}                  # Innovation loading matrix
    iteration_count::S                               # Number of iterations needed for convergence
    node_generator::Function                         # Function to generate the nodes

end

struct SmolyakSolutionDet{T<:Real,S<:Integer} <: ProjectionSolutionDet

    variables::Array{Array{T,1},1}       # Variables
    weights::Array{Array{T,1},1}         # Smolyak weights
    grid::Union{Array{T,2},Array{T,1}}   # Smolyak grid
    multi_index::Array{S,2}              # Smolyak multi index
    layer::Union{S,Array{S,1}}          # Isotropic / anisotropic
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    iteration_count::S                   # Number of iterations needed for convergence
    node_generator::Function             # Function to generate the nodes

end

struct SmolyakSolutionStoch{T<:Real,S<:Integer} <: ProjectionSolutionStoch

    variables::Array{Array{T,1},1}       # Variables
    weights::Array{Array{T,1},1}         # Smolyak weights
    scale_factor::Array{T,1}             # Scale factor for computing scaled weights
    grid::Union{Array{T,2},Array{T,1}}   # Smolyak grid
    multi_index::Array{S,2}              # Smolyak multi index
    layer::Union{S,Array{S,1}}          # Isotropic / anisotropic
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    k::Union{Array{T,2},Array{T,1}}      # Innovation loading matrix
    iteration_count::S                   # Number of iterations needed for convergence
    node_generator::Function             # Function to generate the nodes

end

struct HyperbolicCrossSolutionDet{T<:Real,S<:Integer} <: ProjectionSolutionDet

    variables::Array{Array{T,1},1}       # Variables
    weights::Array{Array{T,1},1}         # Hyperbolic cross weights
    grid::Union{Array{T,2},Array{T,1}}   # Hyperbolic cross grid
    multi_index::Array{S,2}              # Hyperbolic cross multi index
    layer::S                             # Number of approximation layers
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    iteration_count::S                   # Number of iterations needed for convergence
    node_generator::Function             # Function to generate the nodes

end

struct HyperbolicCrossSolutionStoch{T<:Real,S<:Integer} <: ProjectionSolutionStoch

    variables::Array{Array{T,1},1}       # Variables
    weights::Array{Array{T,1},1}         # Hyperbolic cross weights
    scale_factor::Array{T,1}             # Scale factor for computing scaled weights
    grid::Union{Array{T,2},Array{T,1}}   # Hyperbolic cross grid
    multi_index::Array{S,2}              # Hyperbolic cross multi index
    layer::S                             # Number of approximation layers
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    k::Union{Array{T,2},Array{T,1}}      # Innovation loading matrix
    iteration_count::S                   # Number of iterations needed for convergence
    node_generator::Function             # Function to generate the nodes

end

struct PiecewiseLinearSolutionDet{T<:Real,S<:Integer,N} <: ProjectionSolutionDet

    variables::Array{Array{T,N},1}       # Variables
    nodes::Array{Array{T,1},1}           # Nodes
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    iteration_count::S                   # Number of iterations needed for convergence

end

struct PiecewiseLinearSolutionStoch{T<:Real,S<:Integer,N} <: ProjectionSolutionStoch

    variables::Array{Array{T,N},1}       # Variables
    nodes::Array{Array{T,1},1}           # Nodes
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    k::Union{Array{T,2},Array{T,1}}      # Innovation loading matrix
    iteration_count::S                   # Number of iterations needed for convergence

end

##################### Introduce the Equilibrium structure ########################

"""
Constains the solution to a nonlinear model in state-space form, with an additional 
equation for forecasting next-period's jump variable:

x(t+1)      = h(x(t))
  y(t)      = g(x(t))
E_{t}y(t+1) = gh(x(t)) = g(h(x(t)))
"""
struct StateSpaceEqm # Augmented to help compute Euler-equation errors

    g::Function    # Decision rules for jumps
    h::Function    # State transition eqn
    gh::Function   # Forecast eqns for next-period's jump variables

end

"""
Contains the summary information returned from the den_haan_marcet() function.

This summary information: the 1%, 5%, and 10% statistics, which follow a 
Chi-square distribution, along with the number of degrees of freedom, can be 
used to apply the Den Haan and Marcet (1994) statistic for assessing solution 
accuracy.
"""
struct DenHaanMarcetStatistic{T<:Real,S<:Integer}

    one_percent::T
    five_percent::T
    ten_percent::T
    degrees_of_freedom::S

end

################# Struct for prior analysis ########################

struct Prior{N}

    dists::NTuple{N,UnivariateDistribution}

end

function prior(args...)

    prior = Prior(args)

end