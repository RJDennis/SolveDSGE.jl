module SolveDSGE

using LinearAlgebra
using ForwardDiff
using NLsolve
using GaussQuadrature
using ChebyshevApprox
using SmolyakApprox
using PiecewiseLinearApprox
using Random
using ThreadPools
using ThreadsX

include("structures.jl")
include("parser_functions.jl")
include("auxiliary_functions.jl")
include("solution_functions.jl")
include("analysis_functions.jl")

export PerturbationScheme,
       ChebyshevSchemeStoch,
       ChebyshevSchemeDet,
       SmolyakSchemeStoch,
       SmolyakSchemeDet,
       PiecewiseLinearSchemeStoch,
       PiecewiseLinearSchemeDet,
       FirstOrderSolutionStoch,
       FirstOrderSolutionDet,
       SecondOrderSolutionStoch,
       SecondOrderSolutionDet,
       ThirdOrderSolutionStoch,
       ThirdOrderSolutionDet,
       ChebyshevSolutionStoch,
       ChebyshevSolutionDet,
       SmolyakSolutionStoch,
       SmolyakSolutionDet,
       SmolyakSolutionStoch,
       PiecewiseLinearSolutionStoch,
       PiecewiseLinearSolutionDet

export process_model,
       retrieve_processed_model,
       assign_parameters,
       compute_steady_state,
       solve_first_order,
       solve_second_order,
       solve_third_order,
       solve_nonlinear,
       solve_model,
       compute_mean,
       simulate,
       impulses,
       approximate_density,
       approximate_distribution,
       compare_solutions

export chebyshev_nodes,
       chebyshev_extrema,
       chebyshev_extended,
       chebyshev_gauss_lobatto,
       clenshaw_curtis_equidistant,
       chebyshev_evaluate,
       smolyak_evaluate,
       piecewise_linear_evaluate,
       chebyshev_weights,
       smolyak_weights

end
