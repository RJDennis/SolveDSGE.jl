SolveDSGE
=========

SolveDSGE is a Julia package for solving Dynamic Stochastic General Equilibrium (DSGE) models.  The package is aimed at macroeconomists interested in first-order-accurate or second-order-accurate solutions to their general equilibrium models.  SolveDSGE offers a broad array of solution methods that can be applied provided the DSGE model can be expressed in one of several standard dynamic representations.

Installation
------------

You can install SolveDSGE by typing in REPL

Pkg.add("SolveDSGE")

After this you will want to get the latest version from master by typing

Pkg.checkout("SolveDSGE")

First-order-accurate methods
----------------------------

SolveDSGE allows rational expectations equilibria to be computed using a variety of methods.  The appropriate method to use depends on how the (linearized) model is expressed.  The following model forms are accommodated:

- Blanchard and Kahn (1980)
- Klein (2000)
- Sims (2001)
- Binder and Pesaran (1995)

Second-order-accurate methods
-----------------------------

SolveDSGE can also solve DSGE models to second-order-accuracy.  The two solution methods that are included are:

- Gomme and Klein (2011)
- Lombardo and Sutherland (2007)

Optimal policy
------------------

In addition to solving rational expectations models, SolveDSGE also computes optimal policies for a range of model forms.  Specifically, The following policies can be computed for linear-quadratic models

- Discretion
- Commitment
- Quasi-commitment
- Timeless-perspective commitment

Further information
------------------- 

Further information on how to use SolveDSGE is contained in the Package Guide included in the repository.
