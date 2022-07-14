SolveDSGE
=========

SolveDSGE is a Julia package for solving Dynamic Stochastic General Equilibrium (DSGE) models.  The package is aimed at macroeconomists interested in projection-based (primarily) and perturbation-based solutions to their nonlinear general equilibrium models.  SolveDSGE offers a broad array of solution methods that can be applied inter-changeably, allowing one solution to form the initialization for another.  The solution algorithms provided allow a wide-range of model to be solved, including models with correlated shocks, models with shocks that have MA components, and models with occasionally binding constraints.  The package also allows parameters to be assigned values "on-the-fly", allowing it to be combined with simulated method of moments and Bayesian/likelihood-based estimation methods, including local filters and particle filters.  The package has options to solve models using sparse grids, including Smolyak methods and Hyperbolic-cross methods, allowing models with a "large" number of state variables to be solved. 

To solve a model using SolveDSGE two files are required.  The model file contains the specifications, parameters, equations, etc for the model that is to be solved.  The solution file reads and processes the model file, specifies a solution scheme, solves the model, and performs any post-solution analysis.

Although I wish I could claim that the package can solve any model you throw at it, I cannot.  Nevertheless, the solution methods provided are intended to help you solve your nonlinear models, and the wide-range of solution methods provided, combined with homotopy (the package allows a solution using one method to be used as an initialization for another), allow many models to be solved.

If you use this package, then please give it a star on github and acknowledge its use appropriately in your work.

Installation
------------

You can install SolveDSGE by typing in REPL:

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

Further information
-------------------
Version 0.5.x introduces methods intended to help solve models with occasionally binding constraints.  Version 0.4.x introduced multi-threading and represented an important update on Version 0.3, which didn't allow any form of parallization.  Version 0.4.x also removes support for the optimal policy routines that were present in Version 0.2.  At some point in the future, SolveDSGE's approach to solving and analyzing optimal policy is likely to be reintroduced based on the model-file/solution-file framework.

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