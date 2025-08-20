SolveDSGE
=========

SolveDSGE is a Julia package for solving Dynamic Stochastic General Equilibrium (DSGE) models.  The package is aimed at macroeconomists interested in projection-based and perturbation-based solutions to their nonlinear general equilibrium models.  SolveDSGE offers a broad array of solution methods that can be applied inter-changeably, allowing one solution to form the initialization for another.  The solution algorithms provided permit a wide-range of models to be solved, including models with correlated shocks, models with shocks that have MA components, and models with occasionally binding constraints.  The package also lets parameters to be assigned values "on-the-fly", allowing it to be combined with simulated method of moments and Bayesian/likelihood-based estimation methods, including local filters and particle filters.  The package has options to solve models using sparse grids, including Smolyak methods and Hyperbolic-cross methods, allowing models with a "large" number of state variables to be solved. 

To solve a model using SolveDSGE two files are required.  The model file contains the specifications, parameters, equations, etc for the model that is to be solved.  The solution file reads and processes the model file, specifies a solution scheme, solves the model, and performs any post-solution analysis.

Although I wish I could claim that the package can solve any model you throw at it, I cannot.  Nevertheless, the solution methods provided are intended to help you solve your nonlinear models, and the wide-range of solution methods provided, combined with homotopy (the package allows a solution using one method to be used as an initialization for another), allow many models to be solved.

If you use this package, then please give it a star on github and acknowledge its use appropriately in your work.

Installation
------------

You can install SolveDSGE.jl by typing in REPL:

```julia
using Pkg
Pkg.add("SolveDSGE")
```

SolveDSGE is supported by several other packages: ForwardDiff.jl, GaussQuadrature.jl, ChebyshevApprox.jl, SmolyakApprox.jl, HyperbolicCrossApprox.jl, PiecewiseLinearApprox.jl, and NLboxsolve.jl.

Perturbation methods
--------------------

SolveDSGE allows rational expectations equilibria to be computed using first-, second-, third- and fourth-order perturbation methods.

Projection methods
------------------

SolveDSGE can solve DSGE models using projection methods that are based on Chebyshev polynomials (complete and tensor-product polynomials), Smolyak polynomials (isotropic and anisotropic grids), Hyperbolic-cross methods (isotropic and anisotropic grids), and piecewise linear approximation.

Example
-------

For the stochastic growth model, we can represent the model through the model file (model_sgm.txt, say):

#= Example --- Stochastic growth model =#

states:

k, z

end

jumps:

c, ce

end

shocks:

ϵ

end

parameters:

β = 0.99

σ = 2.0

δ = 0.015

α = 0.30

ρ = 0.95

sd = 0.01

end

solvers: Any

equations:

k(+1) = (1.0 - δ)*k + exp(z)*k^α - c

c^(-σ) = β*ce(+1)

ce = c^(-σ)*(1.0 - δ + α*exp(z)*k^(α - 1.0))

z(+1) = ρ*z + sd*ϵ

end

Then the solution file that gets run to solve the model could be:

```julia
using SolveDSGE

filename = "model_sgm.txt"
path = joinpath(@__DIR__,filename)

process_model(path)

# If the model file and the solution file are in the same folder

if !occursin("_processed",path)
    model_processed_path = replace(path,".txt" => "_processed.txt")
end

# Hardcode the path to the processed model file if the model file and the solution file are in different folders

#model_processed_path = "your_path_here\\model_sgm_processed.txt"

include(model_processed_path)

dsge_sgm = retrieve_processed_model()
#dsge_sgm = create_model_structure() # Does the same as the retrieve_processed_model() function

x = [0.0, 34.6, 2.4, 0.2]

tol = 1e-8
maxiters = 1000
ss_obj = compute_steady_state(dsge_sgm,x,tol,maxiters)
ss = ss_obj.zero

P1 = PerturbationScheme(ss,1.0,"first")
P2 = PerturbationScheme(ss,1.0,"second")
P3 = PerturbationScheme(ss,1.0,"third")
P4 = PerturbationScheme(ss,1.0,"fourth")

soln_1o = solve_model(dsge_sgm,P1)
soln_2o = solve_model(dsge_sgm,P2)
soln_3o = solve_model(dsge_sgm,P3)
soln_4o = solve_model(dsge_sgm,P4)

C = ChebyshevSchemeStoch(ss,chebyshev_nodes,[21,21],9,3,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)
#soln_nl_cheb  = solve_model(dsge_sgm,C)
soln_nl_cheb = solve_model(dsge_sgm,soln_1o,C)

S = SmolyakSchemeStoch(ss,chebyshev_gauss_lobatto,9,3,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)
#soln_nl_smol = solve_model(dsge_sgm,S)
soln_nl_smol = solve_model(dsge_sgm,soln_2o,S)

H = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,11,5,9,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)
#soln_nl_hcross = solve_model(dsge_sgm,H)
#soln_nl_hcross = solve_model(dsge_sgm,soln_2o,H)
soln_nl_hcross = solve_model(dsge_sgm,soln_nl_smol,H)

M = PiecewiseLinearSchemeStoch(ss,[21,21],9,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)
soln_nl_pwise = solve_model(dsge_sgm,M)

simulated_data     = simulate(soln_nl_cheb,ss[1:2],100000)
pos_imps, neg_imps = impulses(soln_nl_cheb,50,[1],10000)

ee1, ss1 = euler_errors(dsge_sgm,soln_1o,[0.0960769 26.0; -0.0960769 18.0],1000,123456)
ee2, ss2 = euler_errors(dsge_sgm,soln_2o,[0.0960769 26.0; -0.0960769 18.0],1000,123456)
ee3, ss3 = euler_errors(dsge_sgm,soln_3o,[0.0960769 26.0; -0.0960769 18.0],1000,123456)
ee4, ss4 = euler_errors(dsge_sgm,soln_4o,[0.0960769 26.0; -0.0960769 18.0],1000,123456)

eec, ssc = euler_errors(dsge_sgm,soln_nl_cheb,1000,123456)
ees, sss = euler_errors(dsge_sgm,soln_nl_smol,1000,123456)
eep, ssp = euler_errors(dsge_sgm,soln_nl_hcross,1000,123456)
eeh, ssh = euler_errors(dsge_sgm,soln_nl_pwise,1000,123456)

compare_solutions(soln_1o,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
compare_solutions(soln_2o,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
compare_solutions(soln_3o,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
compare_solutions(soln_4o,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
compare_solutions(soln_nl_smol,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
compare_solutions(soln_nl_hcross,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
compare_solutions(soln_nl_pwise,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
```

Further information
-------------------
Examples of how to use SolveDSGE to solve a model are contained in the examples folder.  At the moment there are four examples: a deterministic model, a stochastic growth model, a New Keynesian sticky price model with the zero lower bound, and a heterogeneous agents model.  The example models are intended to be familiar, easy to follow, and easy to replicate.

Although it no-doubt requires much improvement, there is a User Guide that describes in detail the steps that should be taken to solve a model, and documents the relevant functions, solution schemes, and solution structures.  How to simulate data from a solved model, compute impulse response functions, and approximate PDFs and CDFs is also described.

References
----------

The main research papers upon which the solution methods are based are the following:

Andreasen, M., Fernández-Villaverde, J., and J. Rubio-Ramirez, (2017), "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications", Review of Economic Studies, 0, pp. 1---49.

Binning, A., (2013), "Third-Order Approximation of Dynamic Models Without the Use of Tensors", Norges Bank Working Paper 2013---13.

Dennis R., (2021), "Using a hyperbolic cross to solve non-linear macroeconomic models", CAMA Working Paper No. 93/2021

Gomme, P., and P. Klein, (2011), "Second-Order Approximation of Dynamic Models Without the Use of Tensors", Journal of Economic Dynamics and Control, 35, pp. 604---615.

Judd, K., (1992), "Projection Methods for Solving Aggregate Growth Models", Journal of Economic Theory, 58, pp. 410---452.

Judd, K., Maliar, L., Maliar, S., and R. Valero, (2014), "Smolyak Method for Solving Dynamic Economic Models: Lagrange Interpolation, Anisotropic Grid and Adaptive Domain", Journal of Economic Dynamics and Control, 44, pp. 92---123.

Judd, K., Maliar, L., Maliar, S., and I. Tsener, (2017), "How to Solve Dynamic Stochastic Models Computing Expectations just Once", Quantitative Economics, 8, pp. 851---893.

Klein, P., (2000), "Using the Generalized Schur Form to Solve a Multivariate Linear Rational Expectations Model", Journal of Economic Dynamics and Control, 24, pp. 1405---1423.

Kronmal, R., and M. Tarter, (1968), "The Estimation of Probability Densities and Cumulatives by Fourier Series Methods", Journal of the American Statistical Association, 63, 323, pp. 925--952.

Levintal, O., (2017), "Fifth-Order Perturbation Solution to DSGE models", Journal of Economic Dynamics and Control, 80, pp. 1---16.

Potter, S., (2000), "Nonlinear Impulse Response Functions", Journal of Economic Dynamics and Control, 24, pp. 1425---1446.