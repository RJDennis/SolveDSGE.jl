using SolveDSGE

filename = "model1a.txt"
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

P = ChebyshevSchemeStoch(ss,chebyshev_nodes,[21,21],9,3,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)
PP = ChebyshevSchemeStoch(ss,chebyshev_extrema,[21,21],9,4,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)
PPP = ChebyshevSchemeStoch(ss,chebyshev_extended,[71,71],9,6,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)

soln_nla = solve_model(dsge,P,2)
soln_nlb = solve_model(dsge,soln_fo,PP)
soln_nlc = solve_model(dsge,soln_so,PPP)
soln_nld = solve_model(dsge,soln_to,P,2)
soln_nle = solve_model(dsge,soln_nld,PP)
soln_nlf = solve_model(dsge,soln_nle,PPP)

L = SmolyakSchemeStoch(ss,chebyshev_gauss_lobatto,9,4,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)
LL = SmolyakSchemeStoch(ss,clenshaw_curtis_equidistant,9,4,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)
LLL = SmolyakSchemeStoch(ss,chebyshev_gauss_lobatto,9,4,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)

soln_nlg = solve_model(dsge,L,2)
soln_nlh = solve_model(dsge,soln_to,LL,2)
soln_nli = solve_model(dsge,soln_nlh,LLL)
soln_nlj = solve_model(dsge,soln_nlf,LL)

M = PiecewiseLinearSchemeStoch(ss,[21,21],9,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)
MM = PiecewiseLinearSchemeStoch(ss,[31,31],9,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)

soln_nlk = solve_model(dsge,M,2)
soln_nll = solve_model(dsge,soln_to,MM,2)
soln_nlm = solve_model(dsge,soln_nlf,MM)
soln_nln = solve_model(dsge,soln_nlj,MM)
soln_nlo = solve_model(dsge,soln_nln,MM)

H = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,11,5,9,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)
HH = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,11,6,11,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)
HHH = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,11,12,17,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)

soln_nlz = solve_model(dsge,H)
soln_nly = solve_model(dsge,HH,2)
soln_nlx = solve_model(dsge,soln_to,HHH,2)
soln_nlw = solve_model(dsge,soln_nla,H,2)
soln_nlv = solve_model(dsge,soln_nlh,HHH,2)

simulated_data1 = simulate(soln_nlf,ss[1:2],100000)
pos_imps, neg_imps = impulses(soln_nlf,50,[1],10000)

nodesf, f = approximate_density(simulated_data1[3,:],10,1.900,2.065)
nodesF, F = approximate_distribution(simulated_data1[3,:],10,1.900,2.065)

dyn_eqm_fo  = state_space_eqm(soln_fo)
dyn_eqm_so  = state_space_eqm(soln_so)
dyn_eqm_to  = state_space_eqm(soln_to)
dyn_eqm_foo = state_space_eqm(soln_to)
dyn_eqm_c   = state_space_eqm(soln_nla)
dyn_eqm_s   = state_space_eqm(soln_nlg)
dyn_eqm_p   = state_space_eqm(soln_nlk)
dyn_eqm_h   = state_space_eqm(soln_nlz)

ee1, ss1 = euler_errors(dsge,soln_fo,[0.0960769 26.0; -0.0960769 18.0],1000,123456)
ee2, ss2 = euler_errors(dsge,soln_so,[0.0960769 26.0; -0.0960769 18.0],1000,123456)
ee3, ss3 = euler_errors(dsge,soln_to,[0.0960769 26.0; -0.0960769 18.0],1000,123456)
ee4, ss4 = euler_errors(dsge,soln_foo,[0.0960769 26.0; -0.0960769 18.0],1000,123456)

eec, ssc = euler_errors(dsge,soln_nla,1000,123456)
ees, sss = euler_errors(dsge,soln_nlg,1000,123456)
eep, ssp = euler_errors(dsge,soln_nlk,1000,123456)
eeh, ssh = euler_errors(dsge,soln_nlz,1000,123456)

compare_solutions(soln_fo,soln_nla,soln_nla.domain,100_000,123456)
compare_solutions(soln_so,soln_nla,soln_nla.domain,100_000,123456)
compare_solutions(soln_to,soln_nla,soln_nla.domain,100_000,123456)
compare_solutions(soln_foo,soln_nla,soln_nla.domain,100_000,123456)
compare_solutions(soln_nla,soln_nlg,soln_nla.domain,100_000,123456)
compare_solutions(soln_nla,soln_nlk,soln_nla.domain,100_000,123456)
compare_solutions(soln_nla,soln_nlz,soln_nla.domain,100_000,123456)