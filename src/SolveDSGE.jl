module SolveDSGE

using Random
using LinearAlgebra
using ForwardDiff
using NLsolve
using GaussQuadrature
using ChebyshevApprox
using SmolyakApprox
using PiecewiseLinearApprox
using HyperbolicCrossApprox
using ThreadPools
using ThreadsX
using NLboxsolve
using OrderedCollections

include("structures.jl")
include("parser_functions.jl")
include("auxiliary_functions.jl")
include("solution_functions.jl")
include("analysis_functions.jl")

export PerturbationScheme,
       ChebyshevSchemeStoch,
       ChebyshevSchemeDet,
       ChebyshevSchemeOBCStoch,
       ChebyshevSchemeOBCDet,
       SmolyakSchemeStoch,
       SmolyakSchemeDet,
       SmolyakSchemeOBCStoch,
       SmolyakSchemeOBCDet,
       PiecewiseLinearSchemeStoch,
       PiecewiseLinearSchemeDet,
       PiecewiseLinearSchemeOBCStoch,
       PiecewiseLinearSchemeOBCDet,
       HyperbolicCrossSchemeStoch,
       HyperbolicCrossSchemeDet
       HyperbolicCrossSchemeOBCStoch,
       HyperbolicCrossSchemeOBCDet

export FirstOrderSolutionStoch,
       FirstOrderSolutionDet,
       SecondOrderSolutionStoch,
       SecondOrderSolutionDet,
       ThirdOrderSolutionStoch,
       ThirdOrderSolutionDet,
       FourthOrderSolutionStoch,
       FourthOrderSolutionDet,
       ChebyshevSolutionStoch,
       ChebyshevSolutionDet,
       SmolyakSolutionStoch,
       SmolyakSolutionDet,
       PiecewiseLinearSolutionStoch,
       PiecewiseLinearSolutionDet,
       HyperbolicCrossSolutionStoch,
       HyperbolicCrossSolutionDet,
       StateSpaceEqm

export process_model,
       retrieve_processed_model,
       assign_parameters,
       compute_steady_state,
       solve_first_order,
       solve_second_order,
       solve_third_order,
       solve_model,
       compute_mean,
       simulate,
       impulses,
       approximate_density,
       approximate_distribution,
       compare_solutions,
       state_space_eqm,
       euler_errors,
       den_haan_marcet

export chebyshev_nodes,
       chebyshev_extrema,
       chebyshev_extended,
       chebyshev_gauss_lobatto,
       clenshaw_curtis_equidistant,
       chebyshev_evaluate,
       smolyak_evaluate,
       piecewise_linear_evaluate,
       hyperbolic_cross_evaluate,
       chebyshev_weights,
       smolyak_weights,
       hyperbolic_cross_weights

end
