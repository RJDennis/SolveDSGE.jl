using SolveDSGE

filename = "model_aiyagari.txt"
model_path = joinpath(@__DIR__, filename)

process_model(model_path)

dsge_base = retrieve_processed_model(model_path)

ss = [0.0, 5.5, 1.44, 0.71, 0.0, 0.04, 1.18]

ub = [Inf,Inf,Inf,Inf,Inf,Inf,Inf]
lb = [-Inf,0.0,0.0,0.0,-Inf,0.0,0.0]

kbar = 5.6 # Guess for aggregate capital
dsge = assign_parameters(dsge_base,Dict("kbar"=>kbar))

COBC = ChebyshevSchemeOBCStoch(ss,chebyshev_nodes,[21,21],15,[5,5],[0.4 22.0;-0.4 0.0],lb,ub,1e-8,1e-8,1000,:lm_ar)
soln_nlb = solve_model(dsge,COBC,2)

COBC = ChebyshevSchemeOBCStoch(ss,chebyshev_nodes,[31,31],15,[7,7],[0.4 22.0;-0.4 0.0],lb,ub,1e-8,1e-8,1000,:lm_ar)
soln_nlb = solve_model(dsge,soln_nlb,COBC,2)

# Now compute aggregate capital using the bisection method

ktop = 5.7
kbot = 5.4
iters = 0
while true
    dsge = assign_parameters(dsge_base,Dict("kbar"=>kbar))
    soln_nlb = solve_model(dsge,soln_nlb,COBC,2)
    data = ensemble_simulate(soln_nlb,[0.0,kbar],lb,ub,10_000,100)
    knew = (sum(data)/10_000)[2,end]

    if knew > kbar
        kbot = copy(kbar)
    else
        ktop = copy(kbar)
    end

    kbar = (ktop+kbot)/2
    if abs(ktop-kbot) < 1e-4 || iters > 100
        break
    end

    iters +=1
    println((iters, kbar))

end

data = ensemble_simulate(soln_nlb,[0.0,kbar],lb,ub,100_000,100)
capital = Array{Float64,1}(undef,100_000)
for i in eachindex(capital)
    capital[i] = data[i][2,end]
end

using Plots
histogram(capital,normalize = :true)
