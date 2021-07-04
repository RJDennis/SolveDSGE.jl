using SolveDSGE

filename = "model5a.txt"
path5 = joinpath(@__DIR__,filename)
process_model(path5)

processed_filename = "model5a_processed.txt"
processed_path5 =  joinpath(@__DIR__,processed_filename)

dsge5 = retrieve_processed_model(processed_path5)

x5 = [0.05, 0.05, -0.05, 0.3, 0.3, 0.3, 0.005, 0.005]

tol = 1e-8
maxiters = 1000

ss5 = compute_steady_state(dsge5,x5,tol,maxiters)
ss5 = ss5.zero

N   = PerturbationScheme(ss5,1.0,"first")
NN  = PerturbationScheme(ss5,1.0,"second")
NNN = PerturbationScheme(ss5,1.0,"third")

soln_fo5  = solve_model(dsge5,N)
soln_so5  = solve_model(dsge5,NN)
soln_to5  = solve_model(dsge5,NNN)

P = ChebyshevSchemeStoch(ss5,chebyshev_nodes,[11,11,11,11],9,4,[0.03 0.13 0.01 0.34;-0.03 -0.13 -0.021 0.31],tol,1e-5,maxiters)
PP = ChebyshevSchemeStoch(ss5,chebyshev_nodes,[21,21,21,21],9,4,[0.03 0.13 0.01 0.34;-0.03 -0.13 -0.021 0.31],tol,1e-5,maxiters)
PPP = ChebyshevSchemeStoch(ss5,chebyshev_nodes,[31,31,31,31],9,5,[0.03 0.13 0.01 0.34;-0.03 -0.13 -0.021 0.31],tol,1e-5,maxiters)

soln_nl5a = solve_model(dsge5,P)
soln_nl5b = solve_model(dsge5,soln_to5,PP)
soln_nl5c = solve_model(dsge5,soln_nl5b,PPP)

L = SmolyakSchemeStoch(ss5,chebyshev_gauss_lobatto,9,3,[0.03 0.13 0.01 0.34;-0.03 -0.13 -0.021 0.31],tol,1e-5,maxiters)
LL = SmolyakSchemeStoch(ss5,chebyshev_gauss_lobatto,9,4,[0.03 0.13 0.01 0.34;-0.03 -0.13 -0.021 0.31],tol,1e-5,maxiters)
LLL = SmolyakSchemeStoch(ss5,chebyshev_gauss_lobatto,9,5,[0.03 0.13 0.01 0.34;-0.03 -0.13 -0.021 0.31],tol,1e-5,maxiters)

soln_nl5d = solve_model(dsge5,L)
soln_nl5e = solve_model(dsge5,soln_to5,LL)
soln_nl5f = solve_model(dsge5,soln_nl5e,LLL)

M = PiecewiseLinearSchemeStoch(ss5,[11,11,11,11],9,[0.03 0.13 0.01 0.34;-0.03 -0.13 -0.021 0.31],1e-3,1e-5,maxiters)
MM = PiecewiseLinearSchemeStoch(ss5,[21,21,21,21],9,[0.03 0.13 0.01 0.34;-0.03 -0.13 -0.021 0.31],tol,1e-5,maxiters)

soln_nl5g = solve_model(dsge5,M)
soln_nl5h = solve_model(dsge5,soln_to5,MM)
soln_nl5i = solve_model(dsge5,soln_nl5b,MM)

simulated_data5 = simulate(soln_nl5c,ss5[1:4],100)
pos_imps, neg_imps = impulses(soln_nl5c,50,[1],10000)
