using SolveDSGE

filename = "model1b.txt"
path = joinpath(@__DIR__,filename)
process_model(path)
dsge = retrieve_processed_model(path)

x = [0.05, 21.4, 1.97, 2.8]

tol = 1e-8
maxiters = 1000
ss_obj = compute_steady_state(dsge,x,tol,maxiters)
ss = ss_obj.zero

N    = PerturbationScheme(ss,1.0,"first")
NN   = PerturbationScheme(ss,1.0,"second")
NNN  = PerturbationScheme(ss,1.0,"third")
NNNN = PerturbationScheme(ss,1.0,"fourth")

soln_fo  = solve_model(dsge,N)
soln_so  = solve_model(dsge,NN)
soln_to  = solve_model(dsge,NNN)
soln_foo = solve_model(dsge,NNNN)

P = ChebyshevSchemeDet(ss,chebyshev_nodes,[21,21],3,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)
PP = ChebyshevSchemeDet(ss,chebyshev_nodes,[21,21],4,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)
PPP = ChebyshevSchemeDet(ss,chebyshev_nodes,[71,71],6,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)

soln_nla = solve_model(dsge,P)
soln_nlb = solve_model(dsge,soln_fo,PP)
soln_nlc = solve_model(dsge,soln_so,PPP,2)
soln_nld = solve_model(dsge,soln_to,P)
soln_nle = solve_model(dsge,soln_nld,PP)
soln_nlf = solve_model(dsge,soln_nle,PPP,2)

L = SmolyakSchemeDet(ss,chebyshev_gauss_lobatto,4,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)
LL = SmolyakSchemeDet(ss,clenshaw_curtis_equidistant,4,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)
LLL = SmolyakSchemeDet(ss,chebyshev_gauss_lobatto,5,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)

soln_nlg = solve_model(dsge,L)
soln_nlh = solve_model(dsge,soln_to,LL)
soln_nli = solve_model(dsge,soln_nlh,LLL,2)
soln_nlj = solve_model(dsge,soln_nlf,LL)

M = PiecewiseLinearSchemeDet(ss,[21,21],[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)
MM = PiecewiseLinearSchemeDet(ss,[31,31],[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)

soln_nlk = solve_model(dsge,M)
soln_nll = solve_model(dsge,soln_to,MM)
soln_nlm = solve_model(dsge,soln_nlf,MM,2)
soln_nln = solve_model(dsge,soln_nlj,MM)
soln_nlo = solve_model(dsge,soln_nln,MM)

H = HyperbolicCrossSchemeDet(ss,chebyshev_nodes,5,9,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)
HH = HyperbolicCrossSchemeDet(ss,chebyshev_nodes,6,11,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)
HHH = HyperbolicCrossSchemeDet(ss,chebyshev_nodes,12,17,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)

soln_nlz = solve_model(dsge,H)
soln_nly = solve_model(dsge,HH,2)
soln_nlx = solve_model(dsge,soln_to,HHH,2)
soln_nlw = solve_model(dsge,soln_nla,H,2)
soln_nlv = solve_model(dsge,soln_nlh,HHH,2)

simulated_data1 = simulate(soln_nlf,0.95*ss[1:2],100000)
