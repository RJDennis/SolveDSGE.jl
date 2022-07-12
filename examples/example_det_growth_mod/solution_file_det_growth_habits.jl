using SolveDSGE

filename = "model_det_growth_habits.txt"
path = joinpath(@__DIR__,filename)

process_model(path)

dsge = retrieve_processed_model(path)

x = [0.05, 3.05, 0.7, 0.7, 2.8]

tol = 1e-8
maxiters = 100
ss_obj = compute_steady_state(dsge,x,tol,maxiters)
ss = ss_obj.zero

P   = PerturbationScheme(ss,1.0,"first")
PP  = PerturbationScheme(ss,1.0,"second")
PPP = PerturbationScheme(ss,1.0,"third")

soln_fo = solve_model(dsge,P)
soln_so = solve_model(dsge,PP)
soln_to = solve_model(dsge,PPP)

C   = ChebyshevSchemeDet(ss,chebyshev_nodes,[11,11,11],4,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters,:newton)
CC  = ChebyshevSchemeDet(ss,chebyshev_nodes,[21,21,21],4,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters,:newton)
CCC = ChebyshevSchemeDet(ss,chebyshev_nodes,[31,31,31],5,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters,:newton)

soln_nla = solve_model(dsge,C)
soln_nlb = solve_model(dsge,soln_so,CC)
soln_nlc = solve_model(dsge,soln_nlb,CCC)

S   = SmolyakSchemeDet(ss,chebyshev_gauss_lobatto,3,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters,:newton)
SS  = SmolyakSchemeDet(ss,chebyshev_gauss_lobatto,4,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters,:newton)
SSS = SmolyakSchemeDet(ss,clenshaw_curtis_equidistant,3,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters,:newton)

soln_nld = solve_model(dsge,S)
soln_nle = solve_model(dsge,soln_fo,SS)
soln_nlf = solve_model(dsge,soln_nla,SSS)

H   = HyperbolicCrossSchemeDet(ss,chebyshev_nodes,5,11,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters,:newton)
HH  = HyperbolicCrossSchemeDet(ss,chebyshev_nodes,7,15,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters,:newton)
HHH = HyperbolicCrossSchemeDet(ss,chebyshev_nodes,9,19,[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters,:newton)

soln_nlg = solve_model(dsge,H)
soln_nlh = solve_model(dsge,soln_fo,HH)
soln_nli = solve_model(dsge,soln_nlg,HHH)

M = PiecewiseLinearSchemeDet(ss,[11,11,11],[0.07 3.2 0.8; -0.07 2.9 0.6],tol,1e-6,maxiters,:newton)

soln_nlj = solve_model(dsge,M)
soln_nlk = solve_model(dsge,soln_to,M)
soln_nll = solve_model(dsge,soln_nlc,M)

simulated_data = simulate(soln_nlf,ss[1:3],100_000)
