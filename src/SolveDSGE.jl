module SolveDSGE

using LinearAlgebra
using ForwardDiff
using NLsolve
using GaussQuadrature
using ChebyshevApprox
using SmolyakApprox
using PiecewiseLinearApprox
using Random

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

export get_model,
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
       chebyshev_gauss_lobatto,
       clenshaw_curtis_equidistant,
       chebyshev_evaluate,
       smolyak_evaluate,
       piecewise_linear_evaluate,
       chebyshev_weights,
       smolyak_weights

#### Old functions to retain backward compatability

include("old/types.jl")
#include("old/tracem.jl")
include("old/permutation.jl")
include("old/convert_second_order.jl")
include("old/vec_to_vech.jl")
include("old/vech_to_vec.jl")
include("old/derivative.jl")
include("old/hessian.jl")
include("old/newton.jl")
include("old/doubling.jl")
#include("old/dlyap.jl")
include("old/check_model.jl")
include("old/solve_re.jl")
include("old/solve_disc.jl")
include("old/solve_commit.jl")
include("old/solve_quasi.jl")
include("old/solve_timeless.jl")
include("old/impulses.jl")

export derivative,
       hessian,
       newton,
       solve_re,
       check_model_form,
       convert_second_order

export Blanchard_Kahn_Form,
       Blanchard_Kahn_Soln,
       Klein_Form,
       Klein_Soln,
       Binder_Pesaran_Form,
       Binder_Pesaran_Soln,
       Sims_Form,
       Sims_Soln,
       Gomme_Klein_Form,
       Gomme_Klein_Soln,
       Lombardo_Sutherland_Form,
       Lombardo_Sutherland_Soln

export doubling,
       dlyap,
       solve_disc,
       solve_commit,
       solve_quasi,
       solve_timeless

export State_Space_Objective,
       Structural_Objective

export State_Space_Form,
       Generalized_State_Space_Form,
       Structural_Form,
       Generalized_Structural_Form,
       State_Space_Soln,
       Structural_Soln

export Perturbable_Soln,
       Second_Order_State_Space_Soln

#export impulses

end
