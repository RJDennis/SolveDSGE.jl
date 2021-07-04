using SolveDSGE

filename = "model2b.txt"
path = joinpath(@__DIR__,filename)
process_model(path)
dsge = retrieve_processed_model(path)
dsge = assign_parameters(dsge,[0.99])

x = [0.05, 3.05, 0.7]

tol = 1e-8
maxiters = 1000
ss = compute_steady_state(dsge,x,tol,maxiters)
ss = ss.zero

N   = PerturbationScheme(ss,1.0,"first")
NN  = PerturbationScheme(ss,1.0,"second")
NNN = PerturbationScheme(ss,1.0,"third")

soln_fo = solve_model(dsge,N)
soln_so = solve_model(dsge,NN)
soln_to = solve_model(dsge,NNN)

P = ChebyshevSchemeDet(ss,chebyshev_nodes,[21,21],3,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters)
PP = ChebyshevSchemeDet(ss,chebyshev_nodes,[21,21],4,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters)
PPP = ChebyshevSchemeDet(ss,chebyshev_nodes,[71,71],6,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters)

soln_nla = solve_model(dsge,P)
soln_nlb = solve_model(dsge,soln_fo,PP)
soln_nlc = solve_model(dsge,soln_so,PPP,2)
soln_nld = solve_model(dsge,soln_to,P)
soln_nle = solve_model(dsge,soln_nld,PP)
soln_nlf = solve_model(dsge,soln_nle,PPP,2)

L = SmolyakSchemeDet(ss,chebyshev_gauss_lobatto,4,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters)
LL = SmolyakSchemeDet(ss,clenshaw_curtis_equidistant,4,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters)
LLL = SmolyakSchemeDet(ss,chebyshev_gauss_lobatto,5,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters)

soln_nlg = solve_model(dsge,L)
soln_nlh = solve_model(dsge,soln_to,LL)
soln_nli = solve_model(dsge,soln_nlh,LLL,2)
soln_nlj = solve_model(dsge,soln_nlf,LL)

M = PiecewiseLinearSchemeDet(ss,[21,21],[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters)
MM = PiecewiseLinearSchemeDet(ss,[31,31],[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters)

soln_nlk = solve_model(dsge,M)
soln_nll = solve_model(dsge,soln_to,MM)
soln_nlm = solve_model(dsge,soln_nlf,MM,2)
soln_nln = solve_model(dsge,soln_nlj,MM)
soln_nlo = solve_model(dsge,soln_nln,MM)

simulated_data1 = simulate(soln_nlf,0.95*ss[1:2],100000)
