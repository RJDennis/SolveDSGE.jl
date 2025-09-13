using SolveDSGE

filename = "model2b.txt"
path = joinpath(@__DIR__,filename)
process_model(path)

if !occursin("_processed",path)
    model_processed_path = replace(path,".txt" => "_processed.txt")
end

include(model_processed_path)

dsge = retrieve_processed_model(model)
dsge = assign_parameters(dsge,[0.99])

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

C1 = ChebyshevSchemeDet(ss,chebyshev_nodes,[21,21],3,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:newton)
C2 = ChebyshevSchemeDet(ss,chebyshev_nodes,[21,21],4,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:newton)
C3 = ChebyshevSchemeDet(ss,chebyshev_nodes,[71,71],6,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:newton)

soln_nla = solve_model(dsge,C1)
soln_nlb = solve_model(dsge,soln_1o,C2)
soln_nlc = solve_model(dsge,soln_2o,C3,2)
soln_nld = solve_model(dsge,soln_3o,C1)
soln_nle = solve_model(dsge,soln_nld,C2)
soln_nlf = solve_model(dsge,soln_nle,C3,2)

S1 = SmolyakSchemeDet(ss,chebyshev_gauss_lobatto,4,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:newton)
S2 = SmolyakSchemeDet(ss,clenshaw_curtis_equidistant,4,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:newton)
S3 = SmolyakSchemeDet(ss,chebyshev_gauss_lobatto,5,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:newton)

soln_nlg = solve_model(dsge,S1)
soln_nlh = solve_model(dsge,soln_3o,S2)
soln_nli = solve_model(dsge,soln_nlh,S3,2)
soln_nlj = solve_model(dsge,soln_nlf,S2)

M1 = PiecewiseLinearSchemeDet(ss,[21,21],[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:newton)
M2 = PiecewiseLinearSchemeDet(ss,[31,31],[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:newton)

soln_nlk = solve_model(dsge,M1)
soln_nll = solve_model(dsge,soln_3o,M2)
soln_nlm = solve_model(dsge,soln_nlf,M2,2)
soln_nln = solve_model(dsge,soln_nlj,M2)
soln_nlo = solve_model(dsge,soln_nln,M2)

H1 = HyperbolicCrossSchemeDet(ss,chebyshev_nodes,5,9,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:newton)
H2 = HyperbolicCrossSchemeDet(ss,chebyshev_nodes,6,11,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:newton)
H3 = HyperbolicCrossSchemeDet(ss,chebyshev_nodes,12,17,[0.0960769 26.0; -0.0960769 8.0],tol,1e-6,maxiters,:newton)

soln_nlz = solve_model(dsge,H1)
soln_nly = solve_model(dsge,H2,2)
soln_nlx = solve_model(dsge,soln_3o,H3,2)
soln_nlw = solve_model(dsge,soln_nla,H1,2)
soln_nlv = solve_model(dsge,soln_nlh,H3,2)

simulated_data1 = simulate(soln_nlf,0.95*ss[1:2],100000)
