using SolveDSGE

filename = "model1a.txt"
path = joinpath(@__DIR__,filename)
process_model(path)

if !occursin("_processed",path)
    model_processed_path = replace(path,".txt" => "_processed.txt")
end

include(model_processed_path)

dsge = retrieve_processed_model(model)

x = [0.05, 21.4, 1.97, 2.8]

tol = 1e-8
maxiters = 1000
ss_obj = compute_steady_state(dsge,x,tol,maxiters)
ss = ss_obj.zero

P1 = PerturbationScheme(ss,1.0,"first")
P2 = PerturbationScheme(ss,1.0,"second")
P3 = PerturbationScheme(ss,1.0,"third")
P4 = PerturbationScheme(ss,1.0,"fourth")

soln_1o = solve_model(dsge,P1)
soln_2o = solve_model(dsge,P2)
soln_3o = solve_model(dsge,P3)
soln_4o = solve_model(dsge,P4)

C1 = ChebyshevSchemeStoch(ss,chebyshev_nodes,[21,21],9,3,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)
C2 = ChebyshevSchemeStoch(ss,chebyshev_extrema,[21,21],9,4,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)
C3 = ChebyshevSchemeStoch(ss,chebyshev_extended,[71,71],9,6,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)

soln_nla = solve_model(dsge,C1,2)
soln_nlb = solve_model(dsge,soln_1o,C2)
soln_nlc = solve_model(dsge,soln_2o,C3)
soln_nld = solve_model(dsge,soln_3o,C1,2)
soln_nle = solve_model(dsge,soln_nld,C2)
soln_nlf = solve_model(dsge,soln_nle,C3)

S1 = SmolyakSchemeStoch(ss,chebyshev_gauss_lobatto,9,4,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)
S2 = SmolyakSchemeStoch(ss,clenshaw_curtis_equidistant,9,4,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)
S3 = SmolyakSchemeStoch(ss,chebyshev_gauss_lobatto,9,4,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)

soln_nlg = solve_model(dsge,S1,2)
soln_nlh = solve_model(dsge,soln_3o,S2,2)
soln_nli = solve_model(dsge,soln_nlh,S3)
soln_nlj = solve_model(dsge,soln_nlf,S2)

M1 = PiecewiseLinearSchemeStoch(ss,[21,21],9,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:trust_region)

soln_nlk = solve_model(dsge,M1,2)
soln_nll = solve_model(dsge,soln_3o,M1,2)
soln_nlm = solve_model(dsge,soln_nlf,M1)
soln_nln = solve_model(dsge,soln_nlj,M1)
soln_nlo = solve_model(dsge,soln_nln,M1)

H1 = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,11,5,9,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)
H2 = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,11,6,11,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)
H3 = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,11,12,17,[0.0960769 26.0; -0.0960769 18.0],tol,1e-6,maxiters,:newton)

soln_nlz = solve_model(dsge,H1)
soln_nly = solve_model(dsge,H2,2)
soln_nlx = solve_model(dsge,soln_3o,H3,2)
soln_nlw = solve_model(dsge,soln_nla,H1,2)
soln_nlv = solve_model(dsge,soln_nlh,H3,2)

simulated_data1 = simulate(soln_nlf,ss[1:2],100000)
pos_imps, neg_imps = impulses(soln_nlf,50,[1],10000)

nodesf, f = approximate_density(simulated_data1[3,:],10,1.900,2.065)
nodesF, F = approximate_distribution(simulated_data1[3,:],10,1.900,2.065)

dyn_eqm_1o = state_space_eqm(soln_1o)
dyn_eqm_2o = state_space_eqm(soln_2o)
dyn_eqm_3o = state_space_eqm(soln_3o)
dyn_eqm_4o = state_space_eqm(soln_4o)
dyn_eqm_c  = state_space_eqm(soln_nla)
dyn_eqm_s  = state_space_eqm(soln_nlg)
dyn_eqm_p  = state_space_eqm(soln_nlk)
dyn_eqm_h  = state_space_eqm(soln_nlz)

ee1, ss1 = euler_errors(dsge,soln_1o,[0.0960769 26.0; -0.0960769 18.0],1000,123456)
ee2, ss2 = euler_errors(dsge,soln_2o,[0.0960769 26.0; -0.0960769 18.0],1000,123456)
ee3, ss3 = euler_errors(dsge,soln_3o,[0.0960769 26.0; -0.0960769 18.0],1000,123456)
ee4, ss4 = euler_errors(dsge,soln_4o,[0.0960769 26.0; -0.0960769 18.0],1000,123456)

eec, ssc = euler_errors(dsge,soln_nla,1000,123456)
ees, sss = euler_errors(dsge,soln_nlg,1000,123456)
eep, ssp = euler_errors(dsge,soln_nlk,1000,123456)
eeh, ssh = euler_errors(dsge,soln_nlz,1000,123456)

compare_solutions(soln_1o,soln_nla,soln_nla.domain,100_000,123456)
compare_solutions(soln_2o,soln_nla,soln_nla.domain,100_000,123456)
compare_solutions(soln_3o,soln_nla,soln_nla.domain,100_000,123456)
compare_solutions(soln_4o,soln_nla,soln_nla.domain,100_000,123456)
compare_solutions(soln_nla,soln_nlg,soln_nla.domain,100_000,123456)
compare_solutions(soln_nla,soln_nlk,soln_nla.domain,100_000,123456)
compare_solutions(soln_nla,soln_nlz,soln_nla.domain,100_000,123456)