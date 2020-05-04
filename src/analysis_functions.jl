function compute_mean(soln::R) where {R <: PerturbationSolution}

    if typeof(soln) <: FirstOrderSolutionStoch

        mean_states = zeros(length(soln.hbar))
        mean_jumps  = zeros(length(soln.gbar))

    elseif typeof(soln) <: SecondOrderSolutionStoch

        nx = length(soln.hbar)
        ny = length(soln.gbar)
        hxx = Matrix(reshape(soln.hxx',nx*nx,nx)')
        gxx = Matrix(reshape(soln.gxx',nx*nx,ny)')

        mean_states = (I-soln.hx)\((1/2)*soln.hss + (1/2)*hxx*((I-kron(soln.hx,soln.hx))\(kron(soln.k,soln.k)*vec(soln.sigma))))
        mean_jumps  = soln.gx*mean_states + (1/2)*soln.gss + (1/2)*gxx*((I-kron(soln.hx,soln.hx))\(kron(soln.k,soln.k)*vec(soln.sigma)))

    elseif typeof(soln) <: ThirdOrderSolutionStoch

        mean_states = (I-soln.hx)\((1/2)*soln.hss + (1/2)*soln.hxx*((I-kron(soln.hx,soln.hx))\(kron(soln.k,soln.k)*vec(soln.sigma))))
        mean_jumps  = soln.gx*mean_states + (1/2)*soln.gss + (1/2)*soln.gxx*((I-kron(soln.hx,soln.hx))\(kron(soln.k,soln.k)*vec(soln.sigma)))

        term1 = (I-kron(kron(soln.hx,soln.hx),soln.hx))\(kron(kron(soln.k,soln.k),soln.k)*vec(soln.skewness))
        term2 = (I-kron(soln.hx,soln.hx))\(kron(soln.hx,(1/2)*soln.hxx)*term1)

        mean_states += (I-soln.hx)\(soln.hxx*term2 + (1/6)*soln.hxxx*term1 + (1/6)*soln.hsss)
        mean_jumps += soln.gx*((I-soln.hx)\(soln.hxx*term2 + (1/6)*soln.hxxx*term1 + (1/6)*soln.hsss)) + soln.gxx*term2 + (1/6)*soln.gxxx*term1 + (1/6)*soln.gsss

    else # All deterministic cases are handled here

        mean_states = zeros(length(soln.hbar))
        mean_jumps  = zeros(length(soln.gbar))

    end

    return mean_states.+soln.hbar, mean_jumps.+soln.gbar

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R <: FirstOrderSolutionDet, T <: AbstractFloat, S <: Integer}

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    simulated_states_f      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_f       = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] = initial_state - soln.hbar

    for i = 2:sim_length+1
        simulated_states_f[:,i]  = soln.hx*simulated_states_f[:,i-1]
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
    end

    return [simulated_states_f[:,1:sim_length].+soln.hbar;simulated_jumps_f[:,1:end].+soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S;rndseed=123456) where {R <: FirstOrderSolutionStoch, T <: AbstractFloat, S <: Integer}

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    simulated_states_f      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_f       = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] = initial_state - soln.hbar

    for i = 2:sim_length+1
        simulated_states_f[:,i]  = soln.hx*simulated_states_f[:,i-1] + soln.k*randn(size(soln.sigma,1))
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
    end

    return [simulated_states_f[:,1:sim_length].+soln.hbar;simulated_jumps_f[:,1:end].+soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R <: SecondOrderSolutionDet, T <: AbstractFloat, S <: Integer}

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    hxx = Matrix(reshape(soln.hxx',nx*nx,nx)')
    gxx = Matrix(reshape(soln.gxx',nx*nx,ny)')

    simulated_states_f      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_f       = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] = initial_state - soln.hbar

    simulated_states_s      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_s       = Array{T,2}(undef,ny,sim_length)
    simulated_states_s[:,1] = zeros(nx)

    for i = 2:sim_length+1
        simulated_states_f[:,i]  = soln.hx*simulated_states_f[:,i-1]
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
        simulated_states_s[:,i]  = soln.hx*simulated_states_s[:,i-1] + (1/2)*hxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1])
        simulated_jumps_s[:,i-1] = soln.gx*simulated_states_s[:,i-1] + (1/2)*gxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1])
    end

    simulated_states = simulated_states_f + simulated_states_s
    simulated_jumps  = simulated_jumps_f + simulated_jumps_s

    return [simulated_states[:,1:sim_length].+soln.hbar;simulated_jumps[:,1:end].+soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S;rndseed=123456) where {R <: SecondOrderSolutionStoch, T <: AbstractFloat, S <: Integer}

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    hxx = Matrix(reshape(soln.hxx',nx*nx,nx)')
    gxx = Matrix(reshape(soln.gxx',nx*nx,ny)')

    simulated_states_f      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_f       = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] = initial_state - soln.hbar

    simulated_states_s      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_s       = Array{T,2}(undef,ny,sim_length)
    simulated_states_s[:,1] = zeros(nx)

    for i = 2:sim_length+1
        simulated_states_f[:,i]  = soln.hx*simulated_states_f[:,i-1] + soln.k*randn(size(soln.sigma,1))
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
        simulated_states_s[:,i]  = soln.hx*simulated_states_s[:,i-1] + (1/2)*hxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]) + (1/2)*soln.hss
        simulated_jumps_s[:,i-1] = soln.gx*simulated_states_s[:,i-1] + (1/2)*gxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]) + (1/2)*soln.gss
    end

    simulated_states = simulated_states_f + simulated_states_s
    simulated_jumps  = simulated_jumps_f + simulated_jumps_s

    return [simulated_states[:,1:sim_length].+soln.hbar;simulated_jumps[:,1:end].+soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R <: ThirdOrderSolutionDet, T <: AbstractFloat, S <: Integer}

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    simulated_states_f      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_f       = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] = initial_state - soln.hbar

    simulated_states_s      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_s       = Array{T,2}(undef,ny,sim_length)
    simulated_states_s[:,1] = zeros(nx)

    simulated_states_t      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_t       = Array{T,2}(undef,ny,sim_length)
    simulated_states_t[:,1] = zeros(nx)

    for i = 2:sim_length+1
        simulated_states_f[:,i]  = soln.hx*simulated_states_f[:,i-1]
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
        simulated_states_s[:,i]  = soln.hx*simulated_states_s[:,i-1] + (1/2)*soln.hxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1])
        simulated_jumps_s[:,i-1] = soln.gx*simulated_states_s[:,i-1] + (1/2)*soln.gxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1])
        simulated_states_t[:,i]  = soln.hx*simulated_states_t[:,i-1] + (1/2)*soln.hxx*kron(simulated_states_f[:,i-1],simulated_states_s[:,i-1]) + (1/6)*soln.hxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1])
        simulated_jumps_t[:,i-1] = soln.gx*simulated_states_t[:,i-1] + (1/2)*soln.gxx*kron(simulated_states_f[:,i-1],simulated_states_s[:,i-1]) + (1/6)*soln.gxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1])
    end

    simulated_states = simulated_states_f + simulated_states_s + simulated_states_t
    simulated_jumps  = simulated_jumps_f + simulated_jumps_s + simulated_jumps_t

    return [simulated_states[:,1:sim_length].+soln.hbar;simulated_jumps[:,1:end].+soln.gbar]


end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S;rndseed=123456) where {R <: ThirdOrderSolutionStoch, T <: AbstractFloat, S <: Integer}

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    simulated_states_f      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_f       = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] = initial_state - soln.hbar

    simulated_states_s      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_s       = Array{T,2}(undef,ny,sim_length)
    simulated_states_s[:,1] = zeros(nx)

    simulated_states_t      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_t       = Array{T,2}(undef,ny,sim_length)
    simulated_states_t[:,1] = zeros(nx)

    for i = 2:sim_length+1
        simulated_states_f[:,i]  = soln.hx*simulated_states_f[:,i-1] + soln.k*randn(size(soln.sigma,1))
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
        simulated_states_s[:,i]  = soln.hx*simulated_states_s[:,i-1] + (1/2)*soln.hxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]) + (1/2)*soln.hss
        simulated_jumps_s[:,i-1] = soln.gx*simulated_states_s[:,i-1] + (1/2)*soln.gxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]) + (1/2)*soln.gss
        simulated_states_t[:,i]  = soln.hx*simulated_states_t[:,i-1] + (1/2)*soln.hxx*kron(simulated_states_f[:,i-1],simulated_states_s[:,i-1]) + (1/6)*soln.hxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1]) + (3/6)*soln.hssx*simulated_states_f[:,i-1] + (1/6)*soln.hsss
        simulated_jumps_t[:,i-1] = soln.gx*simulated_states_t[:,i-1] + (1/2)*soln.gxx*kron(simulated_states_f[:,i-1],simulated_states_s[:,i-1]) + (1/6)*soln.gxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1]) + (3/6)*soln.gssx*simulated_states_f[:,i-1] + (1/6)*soln.gsss
    end

    simulated_states = simulated_states_f + simulated_states_s + simulated_states_t
    simulated_jumps  = simulated_jumps_f + simulated_jumps_s + simulated_jumps_t

    return [simulated_states[:,1:sim_length].+soln.hbar;simulated_jumps[:,1:end].+soln.gbar]


end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R <: ChebyshevSolutionDet, T <: AbstractFloat, S <: Integer}

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ny = nv - nx

    N = ndims(soln.weights[1])

    w = Array{Array{T,N},1}(undef,length(soln.variables))
    for i = 1:nv
        w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
    end

    simulated_states = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps  = Array{T,2}(undef,ny,sim_length)
    simulated_states[:,1] = initial_state

    for i = 2:sim_length+1
        for j = 1:nx
            simulated_states[j,i] = chebyshev_evaluate(w[j],simulated_states[:,i-1],soln.order,soln.domain)
        end
        for j = 1:ny
            simulated_jumps[j,i-1] = chebyshev_evaluate(w[nx+j],simulated_states[:,i-1],soln.order,soln.domain)
        end
    end

    return [simulated_states[:,1:sim_length];simulated_jumps[:,1:end]]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S;rndseed=123456) where {R <: ChebyshevSolutionStoch, T <: AbstractFloat, S <: Integer}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.sigma,2)
    ny = nv - nx

    chol_decomp = cholesky(soln.sigma)

    N = ndims(soln.weights[1])

    w = Array{Array{T,N},1}(undef,length(soln.variables))
    for i = 1:nv
        w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
    end

    simulated_states = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps  = Array{T,2}(undef,ny,sim_length)
    simulated_states[:,1] = initial_state

    for i = 2:sim_length+1
        for j = 1:nx
            simulated_states[j,i] = chebyshev_evaluate(w[j],simulated_states[:,i-1],soln.order,soln.domain)
        end
        simulated_states[1:ns,i] = chol_decomp.U*randn(ns)
        for j = 1:ny
            simulated_jumps[j,i-1] = chebyshev_evaluate(w[nx+j],simulated_states[:,i-1],soln.order,soln.domain)
        end
    end

    return [simulated_states[:,1:sim_length];simulated_jumps[:,1:end]]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R <: SmolyakSolutionDet, T <: AbstractFloat, S <: Integer}

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ny = nv - nx

    w = Array{Array{T,1},1}(undef,length(soln.variables))
    for i = 1:nv
        w[i] = smolyak_weights(soln.variables[i],soln.grid,soln.multi_index,soln.domain)
    end

    simulated_states = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps  = Array{T,2}(undef,ny,sim_length)
    simulated_states[:,1] = initial_state

    for i = 2:sim_length+1
        for j = 1:nx
            simulated_states[j,i] = smolyak_evaluate(w[j],simulated_states[:,i-1],soln.multi_index,soln.domain)
        end
        for j = 1:ny
            simulated_jumps[j,i-1] = smolyak_evaluate(w[nx+j],simulated_states[:,i-1],soln.multi_index,soln.domain)
        end
    end

    return [simulated_states[:,1:sim_length];simulated_jumps[:,1:end]]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S;rndseed=123456) where {R <: SmolyakSolutionStoch, T <: AbstractFloat, S <: Integer}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.sigma,2)
    ny = nv - nx

    chol_decomp = cholesky(soln.sigma)

    w = Array{Array{T,1},1}(undef,length(soln.variables))
    for i = 1:nv
        w[i] = smolyak_weights(soln.variables[i],soln.grid,soln.multi_index,soln.domain)
    end

    simulated_states = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps  = Array{T,2}(undef,ny,sim_length)
    simulated_states[:,1] = initial_state

    for i = 2:sim_length+1
        for j = 1:nx
            simulated_states[j,i] = smolyak_evaluate(w[j],simulated_states[:,i-1],soln.multi_index,soln.domain)
        end
        simulated_states[1:ns,i] = chol_decomp.U*randn(ns)
        for j = 1:ny
            simulated_jumps[j,i-1] = smolyak_evaluate(w[nx+j],simulated_states[:,i-1],soln.multi_index,soln.domain)
        end
    end

    return [simulated_states[:,1:sim_length];simulated_jumps[:,1:end]]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R <: PiecewiseLinearSolutionDet, T <: AbstractFloat, S <: Integer}

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ny = nv - nx

    simulated_states = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps  = Array{T,2}(undef,ny,sim_length)
    simulated_states[:,1] = initial_state

    for i = 2:sim_length+1
        for j = 1:nx
            simulated_states[j,i] = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states[:,i-1])
        end
        for j = 1:ny
            simulated_jumps[j,i-1] = piecewise_linear_evaluate(soln.variables[nx+j],soln.nodes,simulated_states[:,i-1])
        end
    end

    return [simulated_states[:,1:sim_length];simulated_jumps[:,1:end]]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S;rndseed=123456) where {R <: PiecewiseLinearSolutionStoch, T <: AbstractFloat, S <: Integer}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.sigma,2)
    ny = nv - nx

    chol_decomp = cholesky(soln.sigma)

    simulated_states = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps  = Array{T,2}(undef,ny,sim_length)
    simulated_states[:,1] = initial_state

    for i = 2:sim_length+1
        for j = 1:nx
            simulated_states[j,i] = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states[:,i-1])
        end
        simulated_states[1:ns,i] = chol_decomp.U*randn(ns)
        for j = 1:ny
            simulated_jumps[j,i-1] = piecewise_linear_evaluate(soln.variables[nx+j],soln.nodes,simulated_states[:,i-1])
        end
    end

    return [simulated_states[:,1:sim_length];simulated_jumps[:,1:end]]

end

function impulses(soln::R,n::S,innovation_to_shock::S,reps::S;rndseed=123456) where {R <: FirstOrderSolutionStoch, S <: Integer}

    if innovation_to_shock > size(soln.k,2)
        error("There is no number $innovation_to_shock shock")
    end

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    simulated_states_pos_f      = zeros(nx,n+1)
    simulated_jumps_pos_f       = zeros(ny,n)
    simulated_states_pos_f[:,1] = soln.k[:,innovation_to_shock]

    simulated_states_neg_f      = zeros(nx,n+1)
    simulated_jumps_neg_f       = zeros(ny,n)
    simulated_states_neg_f[:,1] = -soln.k[:,innovation_to_shock]

    for i = 2:n+1
        simulated_states_pos_f[:,i]  = soln.hx*simulated_states_pos_f[:,i-1]
        simulated_jumps_pos_f[:,i-1] = soln.gx*simulated_states_pos_f[:,i-1]
        simulated_states_neg_f[:,i]  = soln.hx*simulated_states_neg_f[:,i-1]
        simulated_jumps_neg_f[:,i-1] = soln.gx*simulated_states_neg_f[:,i-1]
    end

    return [simulated_states_pos_f[:,1:n];simulated_jumps_pos_f], [simulated_states_neg_f[:,1:n];simulated_jumps_neg_f]

end

function impulses(soln::R,n::S,innovation_to_shock::S,reps::S;rndseed=123456) where {R <: SecondOrderSolutionStoch, S <: Integer}

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    hxx = Matrix(reshape(soln.hxx',nx*nx,nx)')
    gxx = Matrix(reshape(soln.gxx',nx*nx,ny)')

    sample = simulate(soln,soln.hss,5*reps+100)
    innovations = randn(size(soln.k,2),n+1)

    impulses_states_pos = zeros(nx,n+1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n+1)
    impulses_jumps_neg  = zeros(ny,n)

    for j = 1:reps
        simulated_states_pos_f = zeros(nx,n+1)
        simulated_jumps_pos_f  = zeros(ny,n)
        simulated_states_pos_s = zeros(nx,n+1)
        simulated_jumps_pos_s  = zeros(ny,n)

        simulated_states_neg_f = zeros(nx,n+1)
        simulated_jumps_neg_f  = zeros(ny,n)
        simulated_states_neg_s = zeros(nx,n+1)
        simulated_jumps_neg_s  = zeros(ny,n)

        simulated_states_base_f = zeros(nx,n+1)
        simulated_jumps_base_f  = zeros(ny,n)
        simulated_states_base_s = zeros(nx,n+1)
        simulated_jumps_base_s  = zeros(ny,n)

        initial_state = sample[1:nx,rand(101:5*reps+100)]
        simulated_states_pos_f[:,1]  = initial_state + soln.k[:,innovation_to_shock]
        simulated_states_neg_f[:,1]  = initial_state - soln.k[:,innovation_to_shock]
        simulated_states_base_f[:,1] = initial_state

        for i = 2:n+1
            simulated_states_pos_f[:,i]  = soln.hx*simulated_states_pos_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_pos_f[:,i-1] = soln.gx*simulated_states_pos_f[:,i-1]
            simulated_states_pos_s[:,i]  = soln.hx*simulated_states_pos_s[:,i-1] + (1/2)*hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + (1/2)*soln.hss
            simulated_jumps_pos_s[:,i-1] = soln.gx*simulated_states_pos_s[:,i-1] + (1/2)*gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + (1/2)*soln.gss

            simulated_states_neg_f[:,i]  = soln.hx*simulated_states_neg_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_neg_f[:,i-1] = soln.gx*simulated_states_neg_f[:,i-1]
            simulated_states_neg_s[:,i]  = soln.hx*simulated_states_neg_s[:,i-1] + (1/2)*hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + (1/2)*soln.hss
            simulated_jumps_neg_s[:,i-1] = soln.gx*simulated_states_neg_s[:,i-1] + (1/2)*gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + (1/2)*soln.gss

            simulated_states_base_f[:,i]  = soln.hx*simulated_states_base_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_base_f[:,i-1] = soln.gx*simulated_states_base_f[:,i-1]
            simulated_states_base_s[:,i]  = soln.hx*simulated_states_base_s[:,i-1] + (1/2)*hxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + (1/2)*soln.hss
            simulated_jumps_base_s[:,i-1] = soln.gx*simulated_states_base_s[:,i-1] + (1/2)*gxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + (1/2)*soln.gss
        end

        impulses_states_pos += (simulated_states_pos_f + simulated_states_pos_s - simulated_states_base_f - simulated_states_base_s)
        impulses_jumps_pos  += (simulated_jumps_pos_f + simulated_jumps_pos_s - simulated_jumps_base_f - simulated_jumps_base_s)
        impulses_states_neg += (simulated_states_neg_f + simulated_states_neg_s - simulated_states_base_f - simulated_states_base_s)
        impulses_jumps_neg  += (simulated_jumps_neg_f + simulated_jumps_neg_s - simulated_jumps_base_f - simulated_jumps_base_s)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n];impulses_jumps_pos], [impulses_states_neg[:,1:n];impulses_jumps_neg]

end

function impulses(soln::R,n::S,innovation_to_shock::S,reps::S;rndseed=123456) where {R <: ThirdOrderSolutionStoch, S <: Integer}

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    sample = simulate(soln,soln.hss,5*reps+100)
    innovations = randn(size(soln.k,2),n+1)

    impulses_states_pos = zeros(nx,n+1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n+1)
    impulses_jumps_neg  = zeros(ny,n)

    for j = 1:reps
        simulated_states_pos_f = zeros(nx,n+1)
        simulated_jumps_pos_f  = zeros(ny,n)
        simulated_states_pos_s = zeros(nx,n+1)
        simulated_jumps_pos_s  = zeros(ny,n)
        simulated_states_pos_t = zeros(nx,n+1)
        simulated_jumps_pos_t  = zeros(ny,n)

        simulated_states_neg_f = zeros(nx,n+1)
        simulated_jumps_neg_f  = zeros(ny,n)
        simulated_states_neg_s = zeros(nx,n+1)
        simulated_jumps_neg_s  = zeros(ny,n)
        simulated_states_neg_t = zeros(nx,n+1)
        simulated_jumps_neg_t  = zeros(ny,n)

        simulated_states_base_f = zeros(nx,n+1)
        simulated_jumps_base_f  = zeros(ny,n)
        simulated_states_base_s = zeros(nx,n+1)
        simulated_jumps_base_s  = zeros(ny,n)
        simulated_states_base_t = zeros(nx,n+1)
        simulated_jumps_base_t  = zeros(ny,n)

        initial_state = sample[1:nx,rand(101:5*reps+100)]
        simulated_states_pos_f[:,1]  = initial_state + soln.k[:,innovation_to_shock]
        simulated_states_neg_f[:,1]  = initial_state - soln.k[:,innovation_to_shock]
        simulated_states_base_f[:,1] = initial_state

        for i = 2:n+1
            simulated_states_pos_f[:,i]  = soln.hx*simulated_states_pos_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_pos_f[:,i-1] = soln.gx*simulated_states_pos_f[:,i-1]
            simulated_states_pos_s[:,i]  = soln.hx*simulated_states_pos_s[:,i-1] + (1/2)*soln.hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + (1/2)*soln.hss
            simulated_jumps_pos_s[:,i-1] = soln.gx*simulated_states_pos_s[:,i-1] + (1/2)*soln.gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + (1/2)*soln.gss
            simulated_states_pos_t[:,i]  = soln.hx*simulated_states_pos_t[:,i-1] + (1/2)*soln.hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_s[:,i-1]) + (1/6)*soln.hxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + (3/6)*soln.hssx*simulated_states_pos_f[:,i-1] + (1/6)*soln.hsss
            simulated_jumps_pos_t[:,i-1] = soln.gx*simulated_states_pos_t[:,i-1] + (1/2)*soln.gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_s[:,i-1]) + (1/6)*soln.gxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + (3/6)*soln.gssx*simulated_states_pos_f[:,i-1] + (1/6)*soln.gsss

            simulated_states_neg_f[:,i]  = soln.hx*simulated_states_neg_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_neg_f[:,i-1] = soln.gx*simulated_states_neg_f[:,i-1]
            simulated_states_neg_s[:,i]  = soln.hx*simulated_states_neg_s[:,i-1] + (1/2)*soln.hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + (1/2)*soln.hss
            simulated_jumps_neg_s[:,i-1] = soln.gx*simulated_states_neg_s[:,i-1] + (1/2)*soln.gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + (1/2)*soln.gss
            simulated_states_neg_t[:,i]  = soln.hx*simulated_states_neg_t[:,i-1] + (1/2)*soln.hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_s[:,i-1]) + (1/6)*soln.hxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + (3/6)*soln.hssx*simulated_states_neg_f[:,i-1] + (1/6)*soln.hsss
            simulated_jumps_neg_t[:,i-1] = soln.gx*simulated_states_neg_t[:,i-1] + (1/2)*soln.gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_s[:,i-1]) + (1/6)*soln.gxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + (3/6)*soln.gssx*simulated_states_neg_f[:,i-1] + (1/6)*soln.gsss

            simulated_states_base_f[:,i]  = soln.hx*simulated_states_base_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_base_f[:,i-1] = soln.gx*simulated_states_base_f[:,i-1]
            simulated_states_base_s[:,i]  = soln.hx*simulated_states_base_s[:,i-1] + (1/2)*soln.hxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + (1/2)*soln.hss
            simulated_jumps_base_s[:,i-1] = soln.gx*simulated_states_base_s[:,i-1] + (1/2)*soln.gxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + (1/2)*soln.gss
            simulated_states_base_t[:,i]  = soln.hx*simulated_states_base_t[:,i-1] + (1/2)*soln.hxx*kron(simulated_states_base_f[:,i-1],simulated_states_pos_s[:,i-1]) + (1/6)*soln.hxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + (3/6)*soln.hssx*simulated_states_base_f[:,i-1] + (1/6)*soln.hsss
            simulated_jumps_base_t[:,i-1] = soln.gx*simulated_states_base_t[:,i-1] + (1/2)*soln.gxx*kron(simulated_states_base_f[:,i-1],simulated_states_pos_s[:,i-1]) + (1/6)*soln.gxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + (3/6)*soln.gssx*simulated_states_base_f[:,i-1] + (1/6)*soln.gsss
        end

        impulses_states_pos += (simulated_states_pos_f + simulated_states_pos_s + simulated_states_pos_t - simulated_states_base_f - simulated_states_base_s - simulated_states_base_t)
        impulses_jumps_pos  += (simulated_jumps_pos_f + simulated_jumps_pos_s + simulated_jumps_pos_t - simulated_jumps_base_f - simulated_jumps_base_s - simulated_jumps_base_t)
        impulses_states_neg += (simulated_states_neg_f + simulated_states_neg_s + simulated_states_neg_t - simulated_states_base_f - simulated_states_base_s - simulated_states_base_t)
        impulses_jumps_neg  += (simulated_jumps_neg_f + simulated_jumps_neg_s + simulated_jumps_neg_t - simulated_jumps_base_f - simulated_jumps_base_s - simulated_jumps_base_t)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n];impulses_jumps_pos], [impulses_states_neg[:,1:n];impulses_jumps_neg]

end

function impulses(soln::R,n::S,innovation_to_shock::S,reps::S;rndseed=123456) where {R <: ChebyshevSolutionStoch, S <: Integer}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.sigma,2)
    ny = nv - nx

    chol_decomp = cholesky(soln.sigma)

    N = ndims(soln.weights[1])

    w = Array{Array{eltype(soln.domain),N},1}(undef,length(soln.variables))
    for i = 1:nv
        w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
    end

    estimated_steady_state = zeros(nx)
    for i = 1:nx
        estimated_steady_state[i] = soln.variables[i][i]
    end

    sample = simulate(soln,estimated_steady_state,5*reps+100)
    innovations = randn(size(soln.sigma,2),n+1)

    impulses_states_pos = zeros(nx,n+1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n+1)
    impulses_jumps_neg  = zeros(ny,n)

    for l = 1:reps
        simulated_states_pos = zeros(nx,n+1)
        simulated_jumps_pos  = zeros(ny,n)

        simulated_states_neg = zeros(nx,n+1)
        simulated_jumps_neg  = zeros(ny,n)

        simulated_states_base = zeros(nx,n+1)
        simulated_jumps_base  = zeros(ny,n)

        initial_state = sample[1:nx,rand(101:5*reps+100)]
        simulated_states_pos[:,1]    = initial_state
        simulated_states_pos[1:ns,1] += chol_decomp.U[:,innovation_to_shock]
        simulated_states_neg[:,1]    = initial_state
        simulated_states_neg[1:ns,1] -= chol_decomp.U[:,innovation_to_shock]
        simulated_states_base[:,1]   = initial_state

        for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = chebyshev_evaluate(w[j],simulated_states_pos[:,i-1],soln.order,soln.domain)
                simulated_states_neg[j,i]  = chebyshev_evaluate(w[j],simulated_states_neg[:,i-1],soln.order,soln.domain)
                simulated_states_base[j,i] = chebyshev_evaluate(w[j],simulated_states_base[:,i-1],soln.order,soln.domain)
            end
            simulated_states_pos[1:ns,i]  += chol_decomp.U*innovations[:,i]
            simulated_states_neg[1:ns,i]  += chol_decomp.U*innovations[:,i]
            simulated_states_base[1:ns,i] += chol_decomp.U*innovations[:,i]
            for j = 1:ny
                simulated_jumps_pos[j,i-1]  = chebyshev_evaluate(w[nx+j],simulated_states_pos[:,i-1],soln.order,soln.domain)
                simulated_jumps_neg[j,i-1]  = chebyshev_evaluate(w[nx+j],simulated_states_neg[:,i-1],soln.order,soln.domain)
                simulated_jumps_base[j,i-1] = chebyshev_evaluate(w[nx+j],simulated_states_base[:,i-1],soln.order,soln.domain)
            end
        end

        impulses_states_pos += (simulated_states_pos - simulated_states_base)
        impulses_jumps_pos  += (simulated_jumps_pos - simulated_jumps_base)
        impulses_states_neg += (simulated_states_neg - simulated_states_base)
        impulses_jumps_neg  += (simulated_jumps_neg - simulated_jumps_base)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n];impulses_jumps_pos], [impulses_states_neg[:,1:n];impulses_jumps_neg]

end

function impulses(soln::R,n::S,innovation_to_shock::S,reps::S;rndseed=123456) where {R <: SmolyakSolutionStoch, S <: Integer}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.sigma,2)
    ny = nv - nx

    chol_decomp = cholesky(soln.sigma)

    w = Array{Array{eltype(soln.domain),1},1}(undef,length(soln.variables))
    for i = 1:nv
        w[i] = smolyak_weights(soln.variables[i],soln.grid,soln.multi_index,soln.domain)
    end

    estimated_steady_state = zeros(nx)
    for i = 1:nx
        estimated_steady_state[i] = soln.variables[i][i]
    end

    sample = simulate(soln,estimated_steady_state,5*reps+100)
    innovations = randn(size(soln.sigma,2),n+1)

    impulses_states_pos = zeros(nx,n+1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n+1)
    impulses_jumps_neg  = zeros(ny,n)

    for l = 1:reps
        simulated_states_pos = zeros(nx,n+1)
        simulated_jumps_pos  = zeros(ny,n)

        simulated_states_neg = zeros(nx,n+1)
        simulated_jumps_neg  = zeros(ny,n)

        simulated_states_base = zeros(nx,n+1)
        simulated_jumps_base  = zeros(ny,n)

        initial_state = sample[1:nx,rand(101:5*reps+100)]
        simulated_states_pos[:,1]    = initial_state
        simulated_states_pos[1:ns,1] += chol_decomp.U[:,innovation_to_shock]
        simulated_states_neg[:,1]    = initial_state
        simulated_states_neg[1:ns,1] -= chol_decomp.U[:,innovation_to_shock]
        simulated_states_base[:,1]   = initial_state

        for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = smolyak_evaluate(w[j],simulated_states_pos[:,i-1],soln.multi_index,soln.domain)
                simulated_states_neg[j,i]  = smolyak_evaluate(w[j],simulated_states_neg[:,i-1],soln.multi_index,soln.domain)
                simulated_states_base[j,i] = smolyak_evaluate(w[j],simulated_states_base[:,i-1],soln.multi_index,soln.domain)
            end
            simulated_states_pos[1:ns,i]  += chol_decomp.U*innovations[:,i]
            simulated_states_neg[1:ns,i]  += chol_decomp.U*innovations[:,i]
            simulated_states_base[1:ns,i] += chol_decomp.U*innovations[:,i]
            for j = 1:ny
                simulated_jumps_pos[j,i-1]  = smolyak_evaluate(w[nx+j],simulated_states_pos[:,i-1],soln.multi_index,soln.domain)
                simulated_jumps_neg[j,i-1]  = smolyak_evaluate(w[nx+j],simulated_states_neg[:,i-1],soln.multi_index,soln.domain)
                simulated_jumps_base[j,i-1] = smolyak_evaluate(w[nx+j],simulated_states_base[:,i-1],soln.multi_index,soln.domain)
            end
        end

        impulses_states_pos += (simulated_states_pos - simulated_states_base)
        impulses_jumps_pos  += (simulated_jumps_pos - simulated_jumps_base)
        impulses_states_neg += (simulated_states_neg - simulated_states_base)
        impulses_jumps_neg  += (simulated_jumps_neg - simulated_jumps_base)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n];impulses_jumps_pos], [impulses_states_neg[:,1:n];impulses_jumps_neg]

end

function impulses(soln::R,n::S,innovation_to_shock::S,reps::S;rndseed=123456) where {R <: PiecewiseLinearSolutionStoch, S <: Integer}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.sigma,2)
    ny = nv - nx

    chol_decomp = cholesky(soln.sigma)

    estimated_steady_state = vec((soln.domain[1,:] + soln.domain[2,:]))/2

    sample = simulate(soln,estimated_steady_state,5*reps+100)
    innovations = randn(size(soln.sigma,2),n+1)

    impulses_states_pos = zeros(nx,n+1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n+1)
    impulses_jumps_neg  = zeros(ny,n)

    for l = 1:reps
        simulated_states_pos = zeros(nx,n+1)
        simulated_jumps_pos  = zeros(ny,n)

        simulated_states_neg = zeros(nx,n+1)
        simulated_jumps_neg  = zeros(ny,n)

        simulated_states_base = zeros(nx,n+1)
        simulated_jumps_base  = zeros(ny,n)

        initial_state = sample[1:nx,rand(101:5*reps+100)]
        simulated_states_pos[:,1]    = initial_state
        simulated_states_pos[1:ns,1] += chol_decomp.U[:,innovation_to_shock]
        simulated_states_neg[:,1]    = initial_state
        simulated_states_neg[1:ns,1] -= chol_decomp.U[:,innovation_to_shock]
        simulated_states_base[:,1]   = initial_state

        for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_pos[:,i-1])
                simulated_states_neg[j,i]  = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_neg[:,i-1])
                simulated_states_base[j,i] = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_base[:,i-1])
            end
            simulated_states_pos[1:ns,i]  += chol_decomp.U*innovations[:,i]
            simulated_states_neg[1:ns,i]  += chol_decomp.U*innovations[:,i]
            simulated_states_base[1:ns,i] += chol_decomp.U*innovations[:,i]
            for j = 1:ny
                simulated_jumps_pos[j,i-1]  = piecewise_linear_evaluate(soln.variables[nx+j],soln.nodes,simulated_states_pos[:,i-1])
                simulated_jumps_neg[j,i-1]  = piecewise_linear_evaluate(soln.variables[nx+j],soln.nodes,simulated_states_neg[:,i-1])
                simulated_jumps_base[j,i-1] = piecewise_linear_evaluate(soln.variables[nx+j],soln.nodes,simulated_states_base[:,i-1])
            end
        end

        impulses_states_pos += (simulated_states_pos - simulated_states_base)
        impulses_jumps_pos  += (simulated_jumps_pos - simulated_jumps_base)
        impulses_states_neg += (simulated_states_neg - simulated_states_base)
        impulses_jumps_neg  += (simulated_jumps_neg - simulated_jumps_base)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n];impulses_jumps_pos], [impulses_states_neg[:,1:n];impulses_jumps_neg]

end

function approximate_density(sample::Array{T,1},point::T,order::S,a::T,b::T) where {T<:AbstractFloat, S<:Integer}

    n = 0

    c = zeros(order+1)
    for j in eachindex(sample)
        if sample[j] >= a && sample[j] <= b
            n += 1
            for i = 1:order+1
                c[i] += (2.0/(b-a))*cos((i-1)*pi*(sample[j]-a)/(b-a))
            end
        end
    end
    c = c/n

    f = c[1]/2.0
    for i = 2:order+1
        f += c[i]*cos((i-1)*pi*(point-a)/(b-a))
    end

    return f

end

function approximate_density(sample::Array{T,1},order::S,a::T,b::T) where {T<:AbstractFloat, S<:Integer}

    n = 0

    c = zeros(order+1)
    for j in eachindex(sample)
        if sample[j] >= a && sample[j] <= b
            n += 1
            for i = 1:order+1
                c[i] += (2.0/(b-a))*cos((i-1)*pi*(sample[j]-a)/(b-a))
            end
        end
    end
    c = c/n

    points = range(a,b,length=minimum([round(Int,n/100),100]))
    ff = zeros(length(points))
    for j in eachindex(ff)
        f = c[1]/2.0
        for i = 2:order+1
            f += c[i]*cos((i-1)*pi*(points[j]-a)/(b-a))
        end
        ff[j] = f
    end

    return collect(points), ff

end

function approximate_distribution(sample::Array{T,1},point::T,order::S,a::T,b::T) where {T<:AbstractFloat, S<:Integer}

    n = 0

    c = zeros(order+1)
    for j in eachindex(sample)
        if sample[j] >= a && sample[j] <= b
            n += 1
            for i = 1:order+1
                c[i] += (2.0/(b-a))*cos((i-1)*pi*(sample[j]-a)/(b-a))
            end
        end
    end
    c = c/n

    F = (point-a)/(b-a)
    for i = 2:order+1
        F += (c[i]*(b-a)/(pi*(i-1)))*sin((i-1)*pi*(point-a)/(b-a))
    end

    return F

end

function approximate_density(sample::Array{T,2},point::Array{T,1},order::Array{S,1},a::Array{T,1},b::Array{T,1}) where {T<:AbstractFloat,S<:Integer}

    n_vars = size(sample,2)
    C = Array{Array{T,1}}(undef,n_vars)
    n = zeros(n_vars)

    for i = 1:n_vars
        c = zeros(order[i]+1)
        samp = sample[:,i]
        for j in eachindex(samp)
            if samp[j] >= a[i] && samp[j] <= b[i]
                n[i] += 1
                for k = 1:order[i]+1
                    c[k] += (2.0/(b[i]-a[i]))*cos((k-1)*pi*(samp[j]-a[i])/(b[i]-a[i]))
                end
            end
        end
        C[i] = c/n[i]
        C[i][1] = C[i][1]/2.0
    end

    V = Array{Array{T,1}}(undef,n_vars)
    for i = 1:n_vars
        v = Array{T,1}(undef,order[i]+1)
        for j = 1:order[i]+1
            if j == 1
                v[1] = 1.0
            else
                v[j] = cos((j-1)*pi*(point[i]-a[i])/(b[i]-a[i]))
            end
        end
        V[i] = v
    end

    coefs = C[1]
    vars = V[1]
    for i = 2:n_vars
        coefs = kron(coefs,C[i])
        vars = kron(vars,V[i])
    end

    return (coefs'vars)[1]

end

function approximate_distribution(sample::Array{T,1},order::S,a::T,b::T) where {T<:AbstractFloat, S<:Integer}

    n = 0

    c = zeros(order+1)
    for j in eachindex(sample)
        if sample[j] >= a && sample[j] <= b
            n += 1
            for i = 1:order+1
                c[i] += (2.0/(b-a))*cos((i-1)*pi*(sample[j]-a)/(b-a))
            end
        end
    end
    c = c/n

    points = range(a,b,length=minimum([round(Int,n/100),100]))
    FF = zeros(length(points))
    for j in eachindex(FF)
        F = (points[j]-a)/(b-a)
        for i = 2:order+1
            F += (c[i]*(b-a)/(pi*(i-1)))*sin((i-1)*pi*(points[j]-a)/(b-a))
        end
        FF[j] = F
    end

    return collect(points), FF

end
