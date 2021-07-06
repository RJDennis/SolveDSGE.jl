using SolveDSGE

filename = "model3b.txt"
path3 = joinpath(@__DIR__,filename)
process_model(path3)

processed_filename = "model3b_processed.txt"
processed_path3 =  joinpath(@__DIR__,processed_filename)

dsge3 = retrieve_processed_model(processed_path3)

x3 = [0.05, 3.05, 0.7, 0.7, 2.8]

tol = 1e-8
maxiters = 1000

ss3_obj = compute_steady_state(dsge3,x3,tol,maxiters)
ss3 = ss3_obj.zero

N   = PerturbationScheme(ss3,1.0,"first")
NN  = PerturbationScheme(ss3,1.0,"second")
NNN = PerturbationScheme(ss3,1.0,"third")

soln_fo3 = solve_model(dsge3,N)
soln_so3 = solve_model(dsge3,NN)
soln_to3 = solve_model(dsge3,NNN)

P = ChebyshevSchemeDet(ss3,chebyshev_nodes,[11,11,11],4,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters)
PP = ChebyshevSchemeDet(ss3,chebyshev_nodes,[21,21,21],4,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters)
PPP = ChebyshevSchemeDet(ss3,chebyshev_nodes,[31,31,31],5,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters)

soln_nl3a = solve_model(dsge3,P)
soln_nl3b = solve_model(dsge3,soln_so3,PP)
soln_nl3c = solve_model(dsge3,soln_nl3b,PPP)

L = SmolyakSchemeDet(ss3,chebyshev_gauss_lobatto,3,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters)
LL = SmolyakSchemeDet(ss3,chebyshev_gauss_lobatto,3,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters)
LLL = SmolyakSchemeDet(ss3,clenshaw_curtis_equidistant,3,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters)

soln_nl3d = solve_model(dsge3,L)
soln_nl3e = solve_model(dsge3,soln_to3,LL)
soln_nl3f = solve_model(dsge3,soln_nl3a,LLL)

M = PiecewiseLinearSchemeDet(ss3,[11,11,11],[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters)

soln_nl3h = solve_model(dsge3,M)
soln_nl3i = solve_model(dsge3,soln_to3,M)
soln_nl3j = solve_model(dsge3,soln_nl3c,M)

simulated_data3 = simulate(soln_nl3f,ss3[1:3],100)
