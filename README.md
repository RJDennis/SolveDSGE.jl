SolveDSGE
=========

SolveDSGE is a Julia package for solving Dynamic Stochastic General Equilibrium (DSGE) models.  The package is aimed at macroeconomists interested in perturbation-based and projection-based solutions to their nonlinear general equilibrium models.  SolveDSGE offers a broad array of solution methods that can be applied inter-changeably, allowing one solution to form the initialization for another.

To solve a model using SolveDSGE two files are required.  The model file contains the specifications, parameters, equations, etc for the model that is to be solved.  The solution file reads and processes the model file, specifies a solution scheme, solves the model, and performs any post-solution analysis.

Installation
------------

You can install SolveDSGE by typing in REPL

```
using Pkg
Pkg.add("SolveDSGE")
```

SolveDSGE is supported by several other packages: NLsolve, ForwardDiff, GaussQuadrature, ChebyshevApprox, SmolyakApprox, and PiecewiseLinearApprox.

Perturbation methods
--------------------

SolveDSGE allows rational expectations equilibria to be computed using first-, second-, and third-order perturbation methods.

Projection methods
------------------

SolveDSGE can also solve DSGE models using projection methods that are based on Chebyshev polynomials, Smolyak polynomials, and piecewise linear approximation.

Further information
-------------------
This version (version 0.3) represents an important update on version 0.2, which didn't allow projection-based solutions, but all of the functionality present in version 0.2 remains.  This includes the optimal policy algorithms.  At some point in the future, SolveDSGE's approach to solving and analyzing optimal policy is likely to be based on the model-file/solution-file framework.  But for now, if you were using version 0.2, your code should still work as before.

Examples of how to use SolveDSGE to solve a model are contained in the examples folder.  At the moment there are three examples: two stochastic models and one deterministic model.  The example models are intended to be familiar, easy to follow, and easy to replicate.

Although it no-doubt requires much improvement, there is a User Guide that describes in detail the steps that should be taken to solve a model, and documents the relevant functions, solution schemes, and solution structures.  How to simulate data from a solved model, compute impulse response functions, and approximate PDFs and CDFs is also described.

References
----------

The main research papers upon which the solution methods are based are the following:

Andreasen, M., Fernández-Villaverde, J., and J. Rubio-Ramirez, (2017), "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications", Review of Economic Studies, 0, pp. 1---49.

Binning, A., (2013), "Third-order approximation of dynamic models without the use of tensors", Norges Bank Working Paper 2013---13.

Gomme, P., and P. Klein, (2011), "Second-Order Approximation of Dynamic Models Without the Use of Tensors", Journal of Economic Dynamics and Control, 35, pp. 604---615.

Judd, K. (1992), "Projection Methods for Solving Aggregate Growth Models", Journal of Economic Theory, 58, pp.410---452.

Judd, K., Maliar, L., Maliar, S., and R. Valero, (2014), "Smolyak Method for Solving Dynamic Economic Models: Lagrange Interpolation, Anisotropic Grid and Adaptive Domain", Journal of Economic Dynamics and Control, 44, pp. 92---123.

Judd, K., Maliar, L., Maliar, S., and I. Tsener, (2017), "How to Solve Dynamic Stochastic Models Computing Expectations just Once", Quantitative Economics, 8, pp.851---893.

Klein, P., (2000), "Using the Generalized Schur Form to Solve a Multivariate Linear Rational Expectations Model", Journal of Economic Dynamics and Control, 24, pp. 1405---1423.

Kronmal, R., and M. Tarter, (1968), "The Estimation of Probability Densities and Cumulatives by Fourier Series Methods", Journal of the American Statistical Association, 63, 323, pp.925--952.

Levintal, O., (2017), "Fifth-Order Perturbation Solution to DSGE models", Journal of Economic Dynamics and Control, 80, pp. 1---16.

Potter, S., (2000), "Nonlinear Impulse Response Functions", Journal of Economic Dynamics and Control, 24, pp. 1425---1446.
