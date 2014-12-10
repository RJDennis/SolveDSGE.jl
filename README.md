SolveDSGE
=========

SolveDSGE is a Julia package for solving Dynamic Stochastic General Equilibrium (DSGE) models.  The package is aimed at macroeconomists interested in first-order-accurate or second-order-accurate solutions to their general equilibrium models.  SolveDSGE offers a broad array of solution methods that can be applied provided the DSGE model can be expressed in one of several standard dynamic representations.

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

Futher information
------------------

Further information on how to use SolveDSGE is contained in the Package Guide included in the respository.
