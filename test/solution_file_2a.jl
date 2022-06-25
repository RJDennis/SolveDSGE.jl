using SolveDSGE

filename = "model2a.txt"
path = joinpath(@__DIR__,filename)
process_model(path)
dsge = retrieve_processed_model(path)
dsge = assign_parameters(dsge,[0.99])

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

P = ChebyshevSchemeStoch(ss,chebyshev_nodes,[21,21],9,3,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)
PP = ChebyshevSchemeStoch(ss,chebyshev_extrema,[21,21],9,4,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)
PPP = ChebyshevSchemeStoch(ss,chebyshev_extended,[71,71],9,6,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)

soln_nla = solve_model(dsge,P)
soln_nlb = solve_model(dsge,soln_fo,PP)
soln_nlc = solve_model(dsge,soln_so,PPP)
soln_nld = solve_model(dsge,soln_to,P,2)
soln_nle = solve_model(dsge,soln_nld,PP,2)
soln_nlf = solve_model(dsge,soln_nle,PPP)

L = SmolyakSchemeStoch(ss,chebyshev_gauss_lobatto,9,4,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)
LL = SmolyakSchemeStoch(ss,clenshaw_curtis_equidistant,9,4,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)
LLL = SmolyakSchemeStoch(ss,chebyshev_gauss_lobatto,9,5,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)

soln_nlg = solve_model(dsge,L,2)
soln_nlh = solve_model(dsge,soln_to,LL)
soln_nli = solve_model(dsge,soln_nlh,LLL)
soln_nlj = solve_model(dsge,soln_nlf,LL,2)

M = PiecewiseLinearSchemeStoch(ss,[21,21],9,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)
MM = PiecewiseLinearSchemeStoch(ss,[31,31],9,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)

soln_nlk = solve_model(dsge,M)
soln_nll = solve_model(dsge,soln_to,MM)
soln_nlm = solve_model(dsge,soln_nlf,MM)
soln_nln = solve_model(dsge,soln_nlj,MM)
soln_nlo = solve_model(dsge,soln_nln,MM,2)

H = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,11,5,9,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)
HH = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,11,6,11,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)
HHH = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,11,12,17,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:trust_region)

soln_nlz = solve_model(dsge,H)
soln_nly = solve_model(dsge,HH,2)
soln_nlx = solve_model(dsge,soln_to,HHH,2)
soln_nlw = solve_model(dsge,soln_nla,H,2)
soln_nlv = solve_model(dsge,soln_nlh,HHH,2)

simulated_data1 = simulate(soln_nlf,ss[1:2],100000)
pos_imps, neg_imps = impulses(soln_nlf,50,[1],10000)

nodesf, f = approximate_density(simulated_data1[3,:],10,1.95,1.99)
nodesF, F = approximate_distribution(simulated_data1[3,:],10,1.95,1.99)
