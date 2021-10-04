################ Introduce the model structures ###############

abstract type DSGEModel end
abstract type ModelPrimatives end

struct REModelPrimatives{Q <: AbstractString} <: ModelPrimatives

    # This structure contains the critial model information extracted from a
    # model file.

    states::Array{Q,1}
    jumps::Array{Q,1}
    shocks::Array{Q,1}
    variables::Array{Q,1}
    parameters::Array{Q,1}
    parametervalues::Array{Q,1}
    equations::Array{Q,1}
    unassigned_parameters::Array{Q,1}

end

struct REModel{S <: Integer, Q <: AbstractString} <: DSGEModel

    # This structure contains information about a model in a form that the
    # model-solvers can work with.

    number_states::S
    number_jumps::S
    number_shocks::S
    number_variables::S
    number_equations::S
    jumps_approximated::Array{S,1}
    eqns_approximated::Array{S,1}
    variables::Array{Q,1}
    nlsolve_static_function::Function
    static_function::Function
    dynamic_function::Function
    each_eqn_function::Array{Function,1}
    closure_function::Function
    closure_function_piecewise::Function

end

struct REModelPartial{S <: Integer, Q <: AbstractString} <: DSGEModel

    # This structure is similar to REModel, but relates to the case
    # where some parameter values have yet to be assigned.

    number_states::S
    number_jumps::S
    number_shocks::S
    number_variables::S
    number_equations::S
    jumps_approximated::Array{S,1}
    eqns_approximated::Array{S,1}
    variables::Array{Q,1}
    nlsolve_static_function::Function
    static_function::Function
    dynamic_function::Function
    each_eqn_function::Array{Function,1}
    closure_function::Function
    closure_function_piecewise::Function
    unassigned_parameters::Array{Q,1}

end

############### Introduce the solutionscheme structures #######################

abstract type SolutionScheme end
abstract type ProjectionScheme <: SolutionScheme end
abstract type ChebyshevScheme <: ProjectionScheme end
abstract type SmolyakScheme <: ProjectionScheme end
abstract type PiecewiseLinearScheme <: ProjectionScheme end

struct PerturbationScheme{T <: Real, Q <: AbstractString} <: SolutionScheme

    steady_state::Array{T,1}
    cutoff::T
    order::Q

end

struct ChebyshevSchemeStoch{T <: Real, S <: Integer} <: ChebyshevScheme

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    node_number::Union{S,Array{S,1}}
    num_quad_nodes::S
    order::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    tol_fix_point_solver::T
    tol_variables::T
    maxiters::S

end

struct ChebyshevSchemeDet{T <: Real, S <: Integer} <: ChebyshevScheme

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    node_number::Union{S,Array{S,1}}
    order::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    tol_fix_point_solver::T
    tol_variables::T
    maxiters::S

end

struct SmolyakSchemeStoch{T <: Real, S <: Integer} <: SmolyakScheme

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    num_quad_nodes::S
    layer::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    tol_fix_point_solver::T
    tol_variables::T
    maxiters::S

end

struct SmolyakSchemeDet{T <: Real, S <: Integer} <: SmolyakScheme

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    layer::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    tol_fix_point_solver::T
    tol_variables::T
    maxiters::S

end

struct PiecewiseLinearSchemeStoch{T <: Real, S <: Integer} <: PiecewiseLinearScheme

    initial_guess::Union{T,Array{T,1}}
    node_number::Union{S,Array{S,1}}
    num_quad_nodes::S
    domain::Union{Array{T,1},Array{T,2}}
    tol_fix_point_solver::T
    tol_variables::T
    maxiters::S

end

struct PiecewiseLinearSchemeDet{T <: Real, S <: Integer} <: PiecewiseLinearScheme

    initial_guess::Union{T,Array{T,1}}
    node_number::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    tol_fix_point_solver::T
    tol_variables::T
    maxiters::S

end

##################### Introduce the solution structures ########################

abstract type ModelSolution end

abstract type PerturbationSolution <: ModelSolution end
abstract type PerturbationSolutionDet <: PerturbationSolution end
abstract type PerturbationSolutionStoch <: PerturbationSolution end

abstract type ProjectionSolution <: ModelSolution end
abstract type ProjectionSolutionDet <: ProjectionSolution end
abstract type ProjectionSolutionStoch <: ProjectionSolution end

struct FirstOrderSolutionStoch{T <: Real,S <: Integer} <: PerturbationSolutionStoch

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

struct FirstOrderSolutionDet{T <: Real,S <: Integer} <: PerturbationSolutionDet

    # x(t+1) = hx*x(t)
    #   y(t) = gx*x(t)

    hbar::Union{T,Array{T,1}}           # steady state values for predetermined variables
    hx::Union{Array{T,2},Array{T,1}}    # Transition matrix for predetermined variables
    gbar::Union{T,Array{T,1}}           # steady state values for nonpredetermined variables
    gx::Union{Array{T,2},Array{T,1}}    # Decision rule matrix linking nonpredetermined variables to predetermined variables
    grc::S                              # Number of eigenvalues greater than cutoff
    soln_type::String                   # "determinate", "indeterminate", or "explosive"

end

struct SecondOrderSolutionStoch{T <: Real,S <: Integer} <: PerturbationSolutionStoch

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

struct SecondOrderSolutionDet{T <: Real,S <: Integer} <: PerturbationSolutionDet

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

struct ThirdOrderSolutionStoch{T <: Real, S <: Integer} <: PerturbationSolutionStoch

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

struct ThirdOrderSolutionDet{T <: Real, S <: Integer} <: PerturbationSolutionDet

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

struct ChebyshevSolutionStoch{T <: Real, S <: Integer, N} <: ProjectionSolutionStoch

    variables::Array{Array{T,N},1}                   # Variables
    weights::Array{Array{T,N},1}                     # Chebyshev weights
    integrals::Union{Array{T},Array{Array{T,1},1}} # Integrals for computing scaled weights
    nodes::Array{Array{T,1},1}                       # Chebyshev nodes
    order::Union{S,Array{S,1}}                       # Complete polynomial / tensor-product
    domain::Union{Array{T,2},Array{T,1}}             # Domain for state variables / state variable
    k::Union{Array{T,2},Array{T,1}}                  # Innovation loading matrix
    iteration_count::S                               # Number of iterations needed for convergence
    node_generator::Function                         # Function to generate the nodes

end

struct ChebyshevSolutionDet{T <: Real, S <: Integer, N} <: ProjectionSolutionDet

    variables::Array{Array{T,N},1}       # Variables
    weights::Array{Array{T,N},1}         # Chebyshev polynomials
    nodes::Array{Array{T,1},1}           # Chebyshev nodes
    order::Union{S,Array{S,1}}           # Complete polynomial / tensor-product
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    iteration_count::S                   # Number of iterations needed for convergence
    node_generator::Function             # Function to generate the nodes

end

struct SmolyakSolutionStoch{T <: Real, S <: Integer} <: ProjectionSolutionStoch

    variables::Array{Array{T,1},1}       # Variables
    weights::Array{Array{T,1},1}         # Smolyak weights
    scale_factor::Array{T,1}             # Scale factor for computing scaled weights
    grid::Union{Array{T,2},Array{T,1}}   # Smolyak grid
    multi_index::Array{S,2}              # Smolyak multi index
    layer::Union{S,Array{S,1}}           # Isotropic / anisotropic
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    k::Union{Array{T,2},Array{T,1}}      # Innovation loading matrix
    iteration_count::S                   # Number of iterations needed for convergence
    node_generator::Function             # Function to generate the nodes

end

struct SmolyakSolutionDet{T <: Real, S <: Integer} <: ProjectionSolutionDet

    variables::Array{Array{T,1},1}       # Variables
    weights::Array{Array{T,1},1}         # Smolyak polynominals
    grid::Union{Array{T,2},Array{T,1}}   # Smolyak grid
    multi_index::Array{S,2}              # Smolyak multi index
    layer::Union{S,Array{S,1}}           # Isotropic / anisotropic
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    iteration_count::S                   # Number of iterations needed for convergence
    node_generator::Function             # Function to generate the nodes

end

struct PiecewiseLinearSolutionStoch{T <: Real, S <: Integer, N} <: ProjectionSolutionStoch

    variables::Array{Array{T,N},1}       # Variables
    nodes::Array{Array{T,1},1}           # Nodes
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    k::Union{Array{T,2},Array{T,1}}      # Innovation loading matrix
    iteration_count::S                   # Number of iterations needed for convergence

end

struct PiecewiseLinearSolutionDet{T <: Real, S <: Integer, N} <: ProjectionSolutionDet

    variables::Array{Array{T,N},1}       # Variables
    nodes::Array{Array{T,1},1}           # Nodes
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    iteration_count::S                   # Number of iterations needed for convergence

end

##################### Introduce the Equilibrium structure ########################

struct StateSpaceEqm

    # x(t+1)      = h(x(t))
    #   y(t)      = g(x(t))
    # E_{t}y(t+1) = gh(x(t))

    g::Function    # Decision rules for jumps
    h::Function    # State transition eqn
    gh::Function   # Forecast eqns for future jumps

end
