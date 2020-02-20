################ Introduce the model structures ###############

abstract type DSGEModel end

mutable struct REModelPrimatives{Q <: AbstractString,T <: AbstractFloat} <: DSGEModel

    # This structure contains the critial model information extracted from a
    # model file.

    states::Array{Q,1}
    jumps::Array{Q,1}
    shocks::Array{Q,1}
    variables::Array{Q,1}
    lag_variables::Array{Q,1}
    lead_variables::Array{Q,1}
    parameters::Array{Q,1}
    parametervalues::Array{T,1}
    equations::Array{Q,1}

end

mutable struct REModel{S <: Integer, Q <: AbstractString} <: DSGEModel

    # This structure contains information about a model in a form that the
    # model-solvers can work with.

    number_states::S
    number_jumps::S
    number_shocks::S
    number_variables::S
    variables::Array{Q,1}
    number_equations::S
    steady_state_equations::Array{Q,1}
    repackaged_equations::Array{Q,1}
    static_function::Function
    nlsolve_static_function::Function
    dynamic_function::Function
    each_eqn_function::Array{Function,1}
    nonlinear_equations::Array{Q,1}
    closure_function::Function
    closure_function_piecewise::Function
    jumps_approximated::Array{S,1}

end

############### Introduce the solutionscheme structures #######################

abstract type SolutionScheme end
abstract type ProjectionScheme <: SolutionScheme end
abstract type ChebyshevScheme <: ProjectionScheme end
abstract type SmolyakScheme <: ProjectionScheme end
abstract type PiecewiseLinearScheme <: ProjectionScheme end

struct PerturbationScheme{T <: AbstractFloat, Q <: AbstractString} <: SolutionScheme

    steady_state::Array{T,1}
    cutoff::T
    order::Q

end

struct ChebyshevSchemeStoch{T <: AbstractFloat, S <: Integer} <: ChebyshevScheme

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

struct ChebyshevSchemeDet{T <: AbstractFloat, S <: Integer} <: ChebyshevScheme

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    node_number::Union{S,Array{S,1}}
    order::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    tol_fix_point_solver::T
    tol_variables::T
    maxiters::S

end

struct SmolyakSchemeStoch{T <: AbstractFloat, S <: Integer} <: SmolyakScheme

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    num_quad_nodes::S
    layer::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    tol_fix_point_solver::T
    tol_variables::T
    maxiters::S

end

struct SmolyakSchemeDet{T <: AbstractFloat, S <: Integer} <: SmolyakScheme

    initial_guess::Union{T,Array{T,1}}
    node_generator::Function
    layer::Union{S,Array{S,1}}
    domain::Union{Array{T,1},Array{T,2}}
    tol_fix_point_solver::T
    tol_variables::T
    maxiters::S

end

struct PiecewiseLinearSchemeStoch{T <: AbstractFloat, S <: Integer} <: PiecewiseLinearScheme

    initial_guess::Union{T,Array{T,1}}
    node_number::Union{S,Array{S,1}}
    num_quad_nodes::S
    domain::Union{Array{T,1},Array{T,2}}
    tol_fix_point_solver::T
    tol_variables::T
    maxiters::S

end

struct PiecewiseLinearSchemeDet{T <: AbstractFloat, S <: Integer} <: PiecewiseLinearScheme

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

struct FirstOrderSolutionStoch{T <: AbstractFloat,S <: Integer} <: PerturbationSolutionStoch

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

struct FirstOrderSolutionDet{T <: AbstractFloat,S <: Integer} <: PerturbationSolutionDet

    # x(t+1) = hx*x(t)
    #   y(t) = gx*x(t)

    hbar::Union{T,Array{T,1}}           # steady state values for predetermined variables
    hx::Union{Array{T,2},Array{T,1}}    # Transition matrix for predetermined variables
    gbar::Union{T,Array{T,1}}           # steady state values for nonpredetermined variables
    gx::Union{Array{T,2},Array{T,1}}    # Decision rule matrix linking nonpredetermined variables to predetermined variables
    grc::S                              # Number of eigenvalues greater than cutoff
    soln_type::String                   # "determinate", "indeterminate", or "explosive"

end

struct SecondOrderSolutionStoch{T <: AbstractFloat,S <: Integer} <: PerturbationSolutionStoch

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

struct SecondOrderSolutionDet{T <: AbstractFloat,S <: Integer} <: PerturbationSolutionDet

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

struct ThirdOrderSolutionStoch{T <: AbstractFloat, S <: Integer} <: PerturbationSolutionStoch

    # x(t+1) = hx*x(t) + (1/2)*hss + (1/2)*hxx*[kron(x(t),x(t)]
    #        + (1/6)*hsss + (1/6)*hssx*[x(t)] + (1/6)*hxxx*[kron(x(t),x(t),x(t))]
    #        + k*v(t+1)

    #   y(t) = gx*x(t) + (1/2)*gss + (1/2)*gxx*[kron(x(t),x(t)]
    #        +(1/6)*gsss + (1/6)*gssx*[x(t)] + (1/6)*gxxx*[kron(x(t),x(t),x(t))]

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

struct ThirdOrderSolutionDet{T <: AbstractFloat, S <: Integer} <: PerturbationSolutionDet

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

struct ChebyshevSolutionStoch{T <: AbstractFloat, S <: Integer, N} <: ProjectionSolutionStoch

    variables::Array{Array{T,N},1}       # Variables
    weights::Array{Array{T,N},1}         # Chebyshev polynomials
    nodes::Array{Array{T,1},1}           # Chebyshev nodes
    order::Union{S,Array{S,1}}           # Complete polynomial / tensor-product
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    sigma::Union{Array{T,2},Array{T,1}}  # Innovation variance-covariance matrix
    iteration_count::S                   # Number of iterations needed for convergence

end

struct ChebyshevSolutionDet{T <: AbstractFloat, S <: Integer, N} <: ProjectionSolutionDet

    variables::Array{Array{T,N},1}       # Variables
    weights::Array{Array{T,N},1}         # Chebyshev polynomials
    nodes::Array{Array{T,1},1}           # Chebyshev nodes
    order::Union{S,Array{S,1}}           # Complete polynomial / tensor-product
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    iteration_count::S                   # Number of iterations needed for convergence

end

struct SmolyakSolutionStoch{T <: AbstractFloat, S <: Integer} <: ProjectionSolutionStoch

    variables::Array{Array{T,1},1}       # Variables
    weights::Array{Array{T,1},1}         # Smolyak polynominals
    grid::Array{T,2}                     # Smolyak grid
    multi_index::Array{S,2}              # Smolyak multi index
    layer::Union{S,Array{S,1}}           # Isotropic / anisotropic
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    sigma::Union{Array{T,2},Array{T,1}}  # Innovation variance-covariance matrix
    iteration_count::S                   # Number of iterations needed for convergence

end

struct SmolyakSolutionDet{T <: AbstractFloat, S <: Integer} <: ProjectionSolutionDet

    variables::Array{Array{T,1},1}       # Variables
    weights::Array{Array{T,1},1}         # Smolyak polynominals
    grid::Array{T,2}                     # Smolyak grid
    multi_index::Array{S,2}              # Smolyak multi index
    layer::Union{S,Array{S,1}}           # Isotropic / anisotropic
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    iteration_count::S                   # Number of iterations needed for convergence

end

struct PiecewiseLinearSolutionStoch{T <: AbstractFloat, S <: Integer, N} <: ProjectionSolutionStoch

    variables::Array{Array{T,N},1}       # Variables
    nodes::Array{Array{T,1},1}           # Nodes
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    sigma::Union{Array{T,2},Array{T,1}}  # Innovation variance-covariance matrix
    iteration_count::S                   # Number of iterations needed for convergence

end

struct PiecewiseLinearSolutionDet{T <: AbstractFloat, S <: Integer, N} <: ProjectionSolutionDet

    variables::Array{Array{T,N},1}       # Variables
    nodes::Array{Array{T,1},1}           # Nodes
    domain::Union{Array{T,2},Array{T,1}} # Domain for state variables / state variable
    iteration_count::S                   # Number of iterations needed for convergence

end
