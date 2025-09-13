using SolveDSGE

filename = "model_nkdsge_zlb.txt"
path = joinpath(@__DIR__, filename)

process_model(path)

if !occursin("_processed",path)
    model_processed_path = replace(path,".txt" => "_processed.txt")
end

include(model_processed_path)

dsge_nk_zlb = retrieve_processed_model()

x = [0.0, 0.0, -0.005012541823544286, 0.3286717154015939, 0.3286717154015939, 3.2180866812572306, 0.003996156184454205, 0.008793969849268423, 0.3286717154015939, 0.00375000000001823, 0.0]

tol = 1e-8
maxiters = 100
ss_obj = compute_steady_state(dsge_nk_zlb,x,tol,maxiters)
ss = ss_obj.zero

lb = [-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0.0,0.0,-Inf,0.0]
ub = [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf]

# Solve the model ignoring the zlb constraint

C = ChebyshevSchemeStoch(ss,chebyshev_nodes,[11,11,11,15],11,6,[0.045 0.2 0.02 0.345; -0.045 -0.21 -0.03 0.31],tol,1e-8,maxiters,:newton)
soln_nla = solve_model(dsge_nk_zlb,C,2)
data = simulate(soln_nla,ss[1:4],100_000)

# Use previous solution as initialization to solve with zlb imposed

COBC = ChebyshevSchemeOBCStoch(ss,chebyshev_nodes,[11,11,11,15],11,6,[0.045 0.2 0.02 0.345; -0.045 -0.21 -0.03 0.31],lb,ub,tol,1e-8,maxiters,:lm_ar)
soln_nlb = solve_model(dsge_nk_zlb,soln_nla,COBC,2)
data_obc = simulate(soln_nlb,ss[1:4],lb,ub,100_000)

using Plots

# Nominal interest rate is 8'th variable in the system

h1 = histogram(data[8,:],normalize = :true,bins=50)
h2 = histogram(data_obc[8,:],normalize = :true,bins=50)
plot(h1,h2)
