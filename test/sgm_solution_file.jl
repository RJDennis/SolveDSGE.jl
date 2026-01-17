using SolveDSGE

filename = "sgm_model.txt"
path = joinpath(@__DIR__,filename)

process_model(path)

# If the model file and the solution file are in the same folder

if !occursin("_processed",path)
    model_processed_path = replace(path,".txt" => "_processed.txt")
end

# Hardcode the path to the processed model file if the model file and the solution file are in different folders

#model_processed_path = "your_path_here\\model_sgm_processed.txt"

include(model_processed_path)

dsge_sgm = retrieve_processed_model(model)
#dsge_sgm = create_model_structure() # Does the same as the retrieve_processed_model() function

x = [0.0, 34.6, 2.4, 0.2]

tol = 1e-8
maxiters = 1000
ss_obj = compute_steady_state(dsge_sgm,x,tol,maxiters)
ss = ss_obj.zero

P1 = PerturbationScheme(ss,1.0,"first")
P2 = PerturbationScheme(ss,1.0,"second")
P3 = PerturbationScheme(ss,1.0,"third")
P4 = PerturbationScheme(ss,1.0,"fourth")

soln_1o = solve_model(dsge_sgm,P1)
soln_2o = solve_model(dsge_sgm,P2)
soln_3o = solve_model(dsge_sgm,P3)
soln_4o = solve_model(dsge_sgm,P4)

C = ChebyshevSchemeStoch(ss,chebyshev_nodes,[21,21],9,[7,7],[0.15 42.5;-0.15 28.0],tol,1e-6,maxiters,:newton)
#soln_nl_cheb  = solve_model(dsge_sgm,C)
soln_nl_cheb = solve_model(dsge_sgm,soln_1o,C)

S = SmolyakSchemeStoch(ss,chebyshev_gauss_lobatto,9,3,[0.15 42.5;-0.15 28.0],tol,1e-6,maxiters,:newton)
#soln_nl_smol = solve_model(dsge_sgm,S)
soln_nl_smol = solve_model(dsge_sgm,soln_2o,S)

H = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,11,5,9,[0.15 42.5;-0.15 28.0],tol,1e-6,maxiters,:newton)
#soln_nl_hcross = solve_model(dsge_sgm,H)
#soln_nl_hcross = solve_model(dsge_sgm,soln_2o,H)
soln_nl_hcross = solve_model(dsge_sgm,soln_nl_smol,H)

M = PiecewiseLinearSchemeStoch(ss,[21,21],9,[0.15 42.5;-0.15 28.0],tol,1e-6,maxiters,:newton)
soln_nl_pwise = solve_model(dsge_sgm,M)

simulated_data     = simulate(soln_4o,ss[1:2],100000)
pos_imps, neg_imps = impulses(soln_nl_cheb,50,[1],10000)

ee1, ss1 = euler_errors(dsge_sgm,soln_1o,[0.15 42.5;-0.15 28.0],1000,123456)
ee2, ss2 = euler_errors(dsge_sgm,soln_2o,[0.15 42.5;-0.15 28.0],1000,123456)
ee3, ss3 = euler_errors(dsge_sgm,soln_3o,[0.15 42.5;-0.15 28.0],1000,123456)
ee4, ss4 = euler_errors(dsge_sgm,soln_4o,[0.15 42.5;-0.15 28.0],1000,123456)

eec, ssc = euler_errors(dsge_sgm,soln_nl_cheb,1000,123456)
ees, sss = euler_errors(dsge_sgm,soln_nl_smol,1000,123456)
eeh, ssh = euler_errors(dsge_sgm,soln_nl_hcross,1000,123456)
eep, ssp = euler_errors(dsge_sgm,soln_nl_pwise,1000,123456)

log10(maximum(abs,ee1))
log10(maximum(abs,ee2))
log10(maximum(abs,ee3))
log10(maximum(abs,ee4))
log10(maximum(abs,eec))
log10(maximum(abs,ees))
log10(maximum(abs,eeh))
log10(maximum(abs,eep))

compare_solutions(soln_1o,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
compare_solutions(soln_2o,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
compare_solutions(soln_3o,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
compare_solutions(soln_4o,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
compare_solutions(soln_nl_smol,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
compare_solutions(soln_nl_hcross,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)
compare_solutions(soln_nl_pwise,soln_nl_cheb,soln_nl_cheb.domain,100_000,123456)