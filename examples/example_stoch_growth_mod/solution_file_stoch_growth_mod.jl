using SolveDSGE

filename = "model_stoch_growth.txt"
path = joinpath(@__DIR__,filename)

process_model(path)

if !occursin("_processed",path)
    model_processed_path = replace(path,".txt" => "_processed.txt")
end

include(model_processed_path)

dsge = retrieve_processed_model(model)

x = [0.0, 21.4, 2.0, 0.25]

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

C   = ChebyshevSchemeStoch(ss,chebyshev_nodes,[21,21],9,[5,5],[0.15 27.0; -0.15 17.0],tol,1e-6,maxiters,:newton)
CC  = ChebyshevSchemeStoch(ss,chebyshev_nodes,[21,21],9,[6,6],[0.15 27.0; -0.15 17.0],tol,1e-6,maxiters,:newton)
CCC = ChebyshevSchemeStoch(ss,chebyshev_nodes,[71,71],9,[7,7],[0.15 27.0; -0.15 17.0],tol,1e-6,maxiters,:newton)

soln_nla = solve_model(dsge,soln_fo,C)
soln_nlb = solve_model(dsge,soln_fo,CC)
soln_nlc = solve_model(dsge,soln_so,CCC)

S   = SmolyakSchemeStoch(ss,chebyshev_gauss_lobatto,9,4,[0.15 27.0; -0.15 17.0],tol,1e-6,maxiters,:newton)
SS  = SmolyakSchemeStoch(ss,clenshaw_curtis_equidistant,9,4,[0.15 27.0; -0.15 17.0],tol,1e-6,maxiters,:newton)
SSS = SmolyakSchemeStoch(ss,chebyshev_gauss_lobatto,9,5,[0.15 27.0; -0.15 17.0],tol,1e-6,maxiters,:newton)

soln_nld = solve_model(dsge,S)
soln_nle = solve_model(dsge,soln_to,SS)
soln_nlf = solve_model(dsge,soln_nla,SSS)

H   = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,9,5,11,[0.15 27.0; -0.15 17.0],tol,1e-6,maxiters,:newton)
HH  = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,9,6,13,[0.15 27.0; -0.15 17.0],tol,1e-6,maxiters,:newton)
HHH = HyperbolicCrossSchemeStoch(ss,chebyshev_nodes,9,7,15,[0.15 27.0; -0.15 17.0],tol,1e-6,maxiters,:newton)

soln_nlg = solve_model(dsge,H)
soln_nlh = solve_model(dsge,soln_to,HH)
soln_nli = solve_model(dsge,soln_nla,HHH)

M = PiecewiseLinearSchemeStoch(ss,[21,21],9,[0.15 27.0; -0.15 17.0],tol,1e-6,maxiters,:newton)
MM = PiecewiseLinearSchemeStoch(ss,[31,31],9,[0.15 27.0; -0.15 17.0],tol,1e-6,maxiters,:newton)

soln_nlj = solve_model(dsge,M)
soln_nlk = solve_model(dsge,soln_to,MM)
soln_nll = solve_model(dsge,soln_nlf,MM)

data = simulate(soln_nla,ss[1:2],100_000)
pos_imps, neg_imps = impulses(soln_nle,50,[1],10_000)

nodesf, f = approximate_density(data[3,:],10,1.70,2.25)
nodesF, F = approximate_distribution(data[3,:],10,1.70,2.25)

ee1, ss1 = euler_errors(dsge,soln_fo,[0.15 27.0; -0.15 17.0],1_000,123456)
ee2, ss2 = euler_errors(dsge,soln_so,[0.15 27.0; -0.15 17.0],1_000,123456)
ee3, ss3 = euler_errors(dsge,soln_to,[0.15 27.0; -0.15 17.0],1_000,123456)
eec, ssc = euler_errors(dsge,soln_nlb,1_000,123456)
ees, sss = euler_errors(dsge,soln_nle,1_000,123456)
eeh, ssh = euler_errors(dsge,soln_nlh,1_000,123456)
eep, ssp = euler_errors(dsge,soln_nlk,1_000,123456)

using Plots

s1 = scatter(ssc[1,:],eec[1,:])
s2 = scatter(ssc[2,:],eec[1,:])

plot(s1,s2)

log10.(maximum(abs.(ee1),dims = 2))
log10.(maximum(abs.(ee2),dims = 2))
log10.(maximum(abs.(ee3),dims = 2))
log10.(maximum(abs.(eec),dims = 2))
log10.(maximum(abs.(ees),dims = 2))
log10.(maximum(abs.(eeh),dims = 2))
log10.(maximum(abs.(eep),dims = 2))

dhm1 = den_haan_marcet(dsge,soln_fo,ss,123456)
dhm2 = den_haan_marcet(dsge,soln_so,ss,123456)
dhm3 = den_haan_marcet(dsge,soln_to,ss,123456)
dhmc = den_haan_marcet(dsge,soln_nlb,ss,123456)
dhms = den_haan_marcet(dsge,soln_nle,ss,123456)
dhmh = den_haan_marcet(dsge,soln_nlh,ss,123456)
dhmp = den_haan_marcet(dsge,soln_nlk,ss,123456)

compare_solutions(soln_fo,soln_nlb,soln_nlb.domain,100_000,123456)
compare_solutions(soln_so,soln_nlb,soln_nlb.domain,100_000,123456)
compare_solutions(soln_to,soln_nlb,soln_nlb.domain,100_000,123456)
compare_solutions(soln_nlb,soln_nle,soln_nlb.domain,100_000,123456)
compare_solutions(soln_nlb,soln_nlh,soln_nlb.domain,100_000,123456)
compare_solutions(soln_nlb,soln_nlk,soln_nlb.domain,100_000,123456)