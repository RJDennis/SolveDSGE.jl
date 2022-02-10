##################### Analysis functions #########################

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

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R <: FirstOrderSolutionDet, T <: Real, S <: Integer}

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    simulated_states_f      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_f       = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] = initial_state - soln.hbar

    for i = 2:sim_length+1
        simulated_states_f[:,i]  = soln.hx*simulated_states_f[:,i-1]
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
    end

    return [simulated_states_f[:,1:sim_length].+soln.hbar;simulated_jumps_f[:,1:end].+soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S;rndseed=123456) where {R <: FirstOrderSolutionStoch, T <: Real, S <: Integer}

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    simulated_states_f      = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_f       = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] = initial_state - soln.hbar

    for i = 2:sim_length+1
        simulated_states_f[:,i]  = soln.hx*simulated_states_f[:,i-1] + soln.k*randn(size(soln.sigma,1))
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
    end

    return [simulated_states_f[:,1:sim_length].+soln.hbar;simulated_jumps_f[:,1:end].+soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R <: SecondOrderSolutionDet, T <: Real, S <: Integer}

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

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

function simulate(soln::R,initial_state::Array{T,1},sim_length::S;rndseed=123456) where {R <: SecondOrderSolutionStoch, T <: Real, S <: Integer}

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

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

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R <: ThirdOrderSolutionDet, T <: Real, S <: Integer}

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

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

function simulate(soln::R,initial_state::Array{T,1},sim_length::S;rndseed=123456) where {R <: ThirdOrderSolutionStoch, T <: Real, S <: Integer}

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

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

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R <: ChebyshevSolutionDet, T <: Real, S <: Integer}

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ny = nv - nx

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    N = ndims(soln.weights[1])

    w = Array{Array{T,N},1}(undef,length(soln.variables))
    if soln.node_generator == chebyshev_nodes
        for i = 1:nv
            w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:nv
            w[i] = chebyshev_weights_extrema(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:nv
            w[i] = chebyshev_weights_extended(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
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

function simulate(soln::R,initial_state::Array{T,1},sim_length::S;rndseed=123456) where {R <: ChebyshevSolutionStoch, T <: Real, S <: Integer}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    N = ndims(soln.weights[1])

    w = Array{Array{T,N},1}(undef,length(soln.variables))
    if soln.node_generator == chebyshev_nodes
        for i = 1:nv
            w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:nv
            w[i] = chebyshev_weights_extrema(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:nv
            w[i] = chebyshev_weights_extended(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    end

    simulated_states = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps  = Array{T,2}(undef,ny,sim_length)
    simulated_states[:,1] = initial_state

    for i = 2:sim_length+1
        for j = 1:nx
            simulated_states[j,i] = chebyshev_evaluate(w[j],simulated_states[:,i-1],soln.order,soln.domain)
        end
        simulated_states[1:ns,i] += soln.k*randn(ns)
        for j = 1:ny
            simulated_jumps[j,i-1] = chebyshev_evaluate(w[nx+j],simulated_states[:,i-1],soln.order,soln.domain)
        end
    end

    return [simulated_states[:,1:sim_length];simulated_jumps[:,1:end]]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R <: SmolyakSolutionDet, T <: Real, S <: Integer}

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ny = nv - nx

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

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

function simulate(soln::R,initial_state::Array{T,1},sim_length::S;rndseed=123456) where {R <: SmolyakSolutionStoch, T <: Real, S <: Integer}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

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
        simulated_states[1:ns,i] += soln.k*randn(ns)
        for j = 1:ny
            simulated_jumps[j,i-1] = smolyak_evaluate(w[nx+j],simulated_states[:,i-1],soln.multi_index,soln.domain)
        end
    end

    return [simulated_states[:,1:sim_length];simulated_jumps[:,1:end]]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R <: PiecewiseLinearSolutionDet, T <: Real, S <: Integer}

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ny = nv - nx

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

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

function simulate(soln::R,initial_state::Array{T,1},sim_length::S;rndseed=123456) where {R <: PiecewiseLinearSolutionStoch, T <: Real, S <: Integer}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    simulated_states = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps  = Array{T,2}(undef,ny,sim_length)
    simulated_states[:,1] = initial_state

    for i = 2:sim_length+1
        for j = 1:nx
            simulated_states[j,i] = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states[:,i-1])
        end
        simulated_states[1:ns,i] += soln.k*randn(ns)
        for j = 1:ny
            simulated_jumps[j,i-1] = piecewise_linear_evaluate(soln.variables[nx+j],soln.nodes,simulated_states[:,i-1])
        end
    end

    return [simulated_states[:,1:sim_length];simulated_jumps[:,1:end]]

end

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S;rndseed=123456) where {R <: FirstOrderSolutionStoch, S <: Integer, T <: Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    simulated_states_pos_f      = zeros(nx,n+1)
    simulated_jumps_pos_f       = zeros(ny,n)
    simulated_states_pos_f[:,1] = soln.k*innovation_vector

    simulated_states_neg_f      = zeros(nx,n+1)
    simulated_jumps_neg_f       = zeros(ny,n)
    simulated_states_neg_f[:,1] = -soln.k*innovation_vector

    for i = 2:n+1
        simulated_states_pos_f[:,i]  = soln.hx*simulated_states_pos_f[:,i-1]
        simulated_jumps_pos_f[:,i-1] = soln.gx*simulated_states_pos_f[:,i-1]
        simulated_states_neg_f[:,i]  = soln.hx*simulated_states_neg_f[:,i-1]
        simulated_jumps_neg_f[:,i-1] = soln.gx*simulated_states_neg_f[:,i-1]
    end

    return [simulated_states_pos_f[:,1:n];simulated_jumps_pos_f], [simulated_states_neg_f[:,1:n];simulated_jumps_neg_f]

end

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S;rndseed=123456) where {R <: SecondOrderSolutionStoch, S <: Integer, T <: Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

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
        simulated_states_pos_f[:,1]  = initial_state + soln.k*innovation_vector
        simulated_states_neg_f[:,1]  = initial_state - soln.k*innovation_vector
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

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S;rndseed=123456) where {R <: SecondOrderSolutionStoch, S <: Integer, T <: Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    hxx = Matrix(reshape(soln.hxx',nx*nx,nx)')
    gxx = Matrix(reshape(soln.gxx',nx*nx,ny)')

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

        simulated_states_pos_f[:,1]  = initial_state + soln.k*innovation_vector
        simulated_states_neg_f[:,1]  = initial_state - soln.k*innovation_vector
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

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S;rndseed=123456) where {R <: ThirdOrderSolutionStoch, S <: Integer, T <: Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

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
        simulated_states_pos_f[:,1]  = initial_state + soln.k*innovation_vector
        simulated_states_neg_f[:,1]  = initial_state - soln.k*innovation_vector
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

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S;rndseed=123456) where {R <: ThirdOrderSolutionStoch, S <: Integer, T <: Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

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

        simulated_states_pos_f[:,1]  = initial_state + soln.k*innovation_vector
        simulated_states_neg_f[:,1]  = initial_state - soln.k*innovation_vector
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

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S;rndseed=123456) where {R <: ChebyshevSolutionStoch, S <: Integer, T <: Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    N = ndims(soln.weights[1])

    w = Array{Array{eltype(soln.domain),N},1}(undef,length(soln.variables))
    if soln.node_generator == chebyshev_nodes
        for i = 1:nv
            w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:nv
            w[i] = chebyshev_weights_extrema(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:nv
            w[i] = chebyshev_weights_extended(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    end

    estimated_steady_state = zeros(nx)
    for i = 1:nx
        estimated_steady_state[i] = soln.variables[i][i]
    end

    sample = simulate(soln,estimated_steady_state,5*reps+100)
    innovations = randn(size(soln.k,2),n+1)

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
        simulated_states_pos[1:ns,1] += soln.k*innovation_vector
        simulated_states_neg[:,1]    = initial_state
        simulated_states_neg[1:ns,1] -= soln.k*innovation_vector
        simulated_states_base[:,1]   = initial_state

        for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = chebyshev_evaluate(w[j],simulated_states_pos[:,i-1],soln.order,soln.domain)
                simulated_states_neg[j,i]  = chebyshev_evaluate(w[j],simulated_states_neg[:,i-1],soln.order,soln.domain)
                simulated_states_base[j,i] = chebyshev_evaluate(w[j],simulated_states_base[:,i-1],soln.order,soln.domain)
            end
            simulated_states_pos[1:ns,i]  += soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  += soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] += soln.k*innovations[:,i]
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

    return [impulses_states_pos[:,1:n];impulses_jumps_pos],[impulses_states_neg[:,1:n];impulses_jumps_neg]

end

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S;rndseed=123456) where {R <: ChebyshevSolutionStoch, S <: Integer, T <: Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    N = ndims(soln.weights[1])

    w = Array{Array{eltype(soln.domain),N},1}(undef,length(soln.variables))
    if soln.node_generator == chebyshev_nodes
        for i = 1:nv
            w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:nv
            w[i] = chebyshev_weights_extrema(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:nv
            w[i] = chebyshev_weights_extended(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    end

    innovations = randn(size(soln.k,2),n+1)

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

        simulated_states_pos[:,1]    = initial_state
        simulated_states_pos[1:ns,1] += soln.k*innovation_vector
        simulated_states_neg[:,1]    = initial_state
        simulated_states_neg[1:ns,1] -= soln.k*innovation_vector
        simulated_states_base[:,1]   = initial_state

        for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = chebyshev_evaluate(w[j],simulated_states_pos[:,i-1],soln.order,soln.domain)
                simulated_states_neg[j,i]  = chebyshev_evaluate(w[j],simulated_states_neg[:,i-1],soln.order,soln.domain)
                simulated_states_base[j,i] = chebyshev_evaluate(w[j],simulated_states_base[:,i-1],soln.order,soln.domain)
            end
            simulated_states_pos[1:ns,i]  += soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  += soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] += soln.k*innovations[:,i]
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

    return [impulses_states_pos[:,1:n];impulses_jumps_pos],[impulses_states_neg[:,1:n];impulses_jumps_neg]

end

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S;rndseed=123456) where {R <: SmolyakSolutionStoch, S <: Integer, T <: Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    w = Array{Array{eltype(soln.domain),1},1}(undef,length(soln.variables))
    for i = 1:nv
        w[i] = smolyak_weights(soln.variables[i],soln.grid,soln.multi_index,soln.domain)
    end

    estimated_steady_state = zeros(nx)
    for i = 1:nx
        estimated_steady_state[i] = soln.variables[i][i]
    end

    sample = simulate(soln,estimated_steady_state,5*reps+100)
    innovations = randn(size(soln.k,2),n+1)

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
        simulated_states_pos[1:ns,1] += soln.k*innovation_vector
        simulated_states_neg[:,1]    = initial_state
        simulated_states_neg[1:ns,1] -= soln.k*innovation_vector
        simulated_states_base[:,1]   = initial_state

        for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = smolyak_evaluate(w[j],simulated_states_pos[:,i-1],soln.multi_index,soln.domain)
                simulated_states_neg[j,i]  = smolyak_evaluate(w[j],simulated_states_neg[:,i-1],soln.multi_index,soln.domain)
                simulated_states_base[j,i] = smolyak_evaluate(w[j],simulated_states_base[:,i-1],soln.multi_index,soln.domain)
            end
            simulated_states_pos[1:ns,i]  += soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  += soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] += soln.k*innovations[:,i]
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

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S;rndseed=123456) where {R <: SmolyakSolutionStoch, S <: Integer, T <: Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    w = Array{Array{eltype(soln.domain),1},1}(undef,length(soln.variables))
    for i = 1:nv
        w[i] = smolyak_weights(soln.variables[i],soln.grid,soln.multi_index,soln.domain)
    end

    innovations = randn(size(soln.k,2),n+1)

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

        simulated_states_pos[:,1]    = initial_state
        simulated_states_pos[1:ns,1] += soln.k*innovation_vector
        simulated_states_neg[:,1]    = initial_state
        simulated_states_neg[1:ns,1] -= soln.k*innovation_vector
        simulated_states_base[:,1]   = initial_state

        for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = smolyak_evaluate(w[j],simulated_states_pos[:,i-1],soln.multi_index,soln.domain)
                simulated_states_neg[j,i]  = smolyak_evaluate(w[j],simulated_states_neg[:,i-1],soln.multi_index,soln.domain)
                simulated_states_base[j,i] = smolyak_evaluate(w[j],simulated_states_base[:,i-1],soln.multi_index,soln.domain)
            end
            simulated_states_pos[1:ns,i]  += soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  += soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] += soln.k*innovations[:,i]
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

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S;rndseed=123456) where {R <: PiecewiseLinearSolutionStoch, S <: Integer, T <: Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    estimated_steady_state = vec((soln.domain[1,:] + soln.domain[2,:]))/2

    sample = simulate(soln,estimated_steady_state,5*reps+100)
    innovations = randn(size(soln.k,2),n+1)

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
        simulated_states_pos[1:ns,1] += soln.k*innovation_vector
        simulated_states_neg[:,1]    = initial_state
        simulated_states_neg[1:ns,1] -= soln.k*innovation_vector
        simulated_states_base[:,1]   = initial_state

        for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_pos[:,i-1])
                simulated_states_neg[j,i]  = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_neg[:,i-1])
                simulated_states_base[j,i] = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_base[:,i-1])
            end
            simulated_states_pos[1:ns,i]  += soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  += soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] += soln.k*innovations[:,i]
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

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S;rndseed=123456) where {R <: PiecewiseLinearSolutionStoch, S <: Integer, T <: Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    innovations = randn(size(soln.k,2),n+1)

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

        simulated_states_pos[:,1]    = initial_state
        simulated_states_pos[1:ns,1] += soln.k*innovation_vector
        simulated_states_neg[:,1]    = initial_state
        simulated_states_neg[1:ns,1] -= soln.k*innovation_vector
        simulated_states_base[:,1]   = initial_state

        for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_pos[:,i-1])
                simulated_states_neg[j,i]  = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_neg[:,i-1])
                simulated_states_base[j,i] = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_base[:,i-1])
            end
            simulated_states_pos[1:ns,i]  += soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  += soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] += soln.k*innovations[:,i]
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

function approximate_density(sample::Array{T,1},point::T,order::S,a::T,b::T) where {T<:Real, S<:Integer}

    if a >= b
        error("'a' must be less than 'b'")
    end

    if point < a || point > b
        error("'point' must be between 'a' and 'b'")
    end

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

function approximate_density(sample::Array{T,1},order::S,a::T,b::T) where {T<:Real, S<:Integer}

    if a >= b
        error("'a' must be less than 'b'")
    end

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

function approximate_distribution(sample::Array{T,1},point::T,order::S,a::T,b::T) where {T<:Real, S<:Integer}

    if a >= b
        error("'a' must be less than 'b'")
    end

    if point < a || point > b
        error("'point' must be between 'a' and 'b'")
    end

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

function approximate_distribution(sample::Array{T,1},order::S,a::T,b::T) where {T<:Real, S<:Integer}

    if a >= b
        error("'a' must be less than 'b'")
    end

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

function compare_solutions(solna::R1,solnb::R2,domain::Array{T,2},seed::S = 123456) where {T <: Real, S <:Integer, R1 <: ModelSolution, R2 <: ModelSolution}

    if typeof(solna) <: PerturbationSolution && typeof(solnb) <: PerturbationSolution
        if length(solna.gbar) != length(solnb.gbar) || length(solna.hbar) != length(solnb.hbar)
            error("The solutions are not of the same model.")
        end
        nx = length(solna.hbar)
        ny = length(solna.gbar)
        sup_errors = zeros(ny)
    elseif typeof(solna) <: PerturbationSolution && typeof(solnb) <: ProjectionSolution
        if length(solna.hbar) != size(solnb.domain,2) || length(solna.gbar) != length(solnb.variables) - size(solnb.domain,2)
            error("The solutions are not of the same model.")
        end
        nx = length(solna.hbar)
        ny = length(solna.gbar)
        sup_errors = zeros(ny)
    elseif typeof(solna) <: ProjectionSolution && typeof(solnb) <: PerturbationSolution
        if length(solnb.hbar) != size(solna.domain,2) || length(solnb.gbar) != length(solna.variables) - size(solna.domain,2)
            error("The solutions are not of the same model.")
        end
        nx = length(solnb.hbar)
        ny = length(solnb.gbar)
        sup_errors = zeros(ny)
    else
        if size(solna.domain) != size(solnb.domain) || length(solna.variables) != length(solnb.variables)
            error("The solutions are not of the same model.")
        end
        nx = size(solna.domain,2)
        ny = length(solna.variables)-nx
        sup_errors = zeros(ny)
    end

    solutions = (solna,solnb)
    Random.seed!(seed)

    n = 100000
    state = domain[2,:] .+ rand(nx,n).*(domain[1,:] - domain[2,:])
    vars = [zeros(ny,n),zeros(ny,n)]
    for j = 1:length(solutions)
        soln = solutions[j]
        for k = 1:ny
            if typeof(soln) <: FirstOrderSolutionDet{T,S}
                for i = 1:n
                    vars[j][k,i] = soln.gbar[k] + (soln.gx[k:k,:]*(state[:,i] - soln.hbar))[1]
                end
            elseif typeof(soln) <: FirstOrderSolutionStoch{T,S}
                for i = 1:n
                    vars[j][k,i] = soln.gbar[k] + (soln.gx[k:k,:]*(state[:,1] - soln.hbar))[1]
                end
            elseif typeof(soln) <: SecondOrderSolutionDet{T,S}
                for i = 1:n
                    vars[j][k,i] = soln.gbar[k] + (soln.gx[k:k,:]*(state[:,i] - soln.hbar))[1]
                    vars[j][k,i] += (1/2)*sum(vec(soln.gxx[(k-1)*nx+1:k*nx,:]).*kron((state[:,i] - soln.hbar),(state[:,i] - soln.hbar)))
                end
            elseif typeof(soln) <: SecondOrderSolutionStoch{T,S}
                for i = 1:n
                    vars[j][k,i] = soln.gbar[k] + (soln.gx[k:k,:]*(state[:,i] - soln.hbar))[1]
                    vars[j][k,i] += (1/2)*soln.gss[k] + (1/2)*sum(vec(soln.gxx[(k-1)*nx+1:k*nx,:]).*kron((state[:,i] - soln.hbar),(state[:,i] - soln.hbar)))
                end
            elseif typeof(soln) <: ThirdOrderSolutionDet{T,S}
                for i = 1:n
                    vars[j][k,i] = soln.gbar[k] + (soln.gx[k:k,:]*(state[:,i] - soln.hbar))[1]
                    vars[j][k,i] += (1/2)*((soln.gxx[k:k,:])*kron((state[:,i] - soln.hbar),(state[:,i] - soln.hbar)))[1] + (1/6)*(soln.gxxx[k:k,:]*kron(kron((state[:,i] - soln.hbar),(state[:,i] - soln.hbar)),(state[:,i] - soln.hbar)))[1]
                end
            elseif typeof(soln) <: ThirdOrderSolutionStoch{T,S}
                for i = 1:n
                    vars[j][k,i] = soln.gbar[k] + (soln.gx[k:k,:]*(state[:,i] - soln.hbar))[1]
                    vars[j][k,i] += (1/2)*soln.gss[k] + (1/2)*((soln.gxx[k:k,:])*kron((state[:,i] - soln.hbar),(state[:,i] - soln.hbar)))[1] + (3/6)*(soln.gssx[k:k,:]*(state[:,i] - soln.hbar))[1] + (1/6)*(soln.gxxx[k:k,:]*kron(kron((state[:,i] - soln.hbar),(state[:,i] - soln.hbar)),(state[:,i] - soln.hbar)))[1]
                end
            elseif typeof(soln) <: Union{ChebyshevSolutionDet{T,S},ChebyshevSolutionStoch{T,S}}
                if soln.node_generator == chebyshev_nodes
                    w = chebyshev_weights(soln.variables[nx+k],soln.nodes,soln.order,soln.domain)
                elseif soln.node_generator == chebyshev_extrema
                    w = chebyshev_weights_extrema(soln.variables[nx+k],soln.nodes,soln.order,soln.domain)
                elseif soln.node_generator == chebyshev_extended
                    w = chebyshev_weights_extended(soln.variables[nx+k],soln.nodes,soln.order,soln.domain)
                end
                for i = 1:n
                    vars[j][k,i] = chebyshev_evaluate(w,state[:,i],soln.order,domain)
                end
            elseif typeof(soln) <: Union{SmolyakSolutionDet{T,S},SmolyakSolutionStoch{T,S}}
                w = smolyak_weights(soln.variables[nx+k],soln.grid,soln.multi_index,soln.domain)
                for i = 1:n
                    vars[j][k,i] = smolyak_evaluate(w,state[:,i],soln.multi_index,domain)
                end
            elseif typeof(soln) <: Union{PiecewiseLinearSolutionDet{T,S},PiecewiseLinearSolutionStoch{T,S}}
                for i = 1:n
                    vars[j][k,i] = piecewise_linear_evaluate(soln.variables[nx+k],soln.nodes,state[:,i])
                end
            end
        end
    end
    for j = 1:ny
        sup_errors[j] = maximum(abs,vars[1][j,:]-vars[2][j,:])
    end

    return sup_errors

end

############################################################################

function decision_rule(soln::R) where {R <: Union{FirstOrderSolutionDet,FirstOrderSolutionStoch}}

    function create_decision_rule(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        
        y = soln.gbar + soln.gx*(state-soln.hbar)

        return y
    
    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R <: SecondOrderSolutionDet}

    function create_decision_rule(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)
        ny = length(soln.gbar)

        y = zeros(ny)
        
        for i = 1:ny
            y[i] = soln.gbar[i] + (soln.gx[i:i,:]*(state-soln.hbar))[1] + (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((state - soln.hbar),(state - soln.hbar)))
        end

        return y
    
    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R <: SecondOrderSolutionStoch}

    function create_decision_rule(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)
        ny = length(soln.gbar)

        y = zeros(ny)
        
        for i = 1:ny
            y[i] = soln.gbar[i] + (soln.gx[i:i,:]*(state-soln.hbar))[1] + (1/2)*soln.gss[i] + (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((state - soln.hbar),(state - soln.hbar)))
        end

        return y
    
    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R <:ThirdOrderSolutionDet}

    function create_decision_rule(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)
        ny = length(soln.gbar)

        y = zeros(ny)
        
        for i = 1:ny
            y[i] = soln.gbar[i] + (soln.gx[i:i,:]*(state-soln.hbar))[1] + (1/2)*((soln.gxx[i:i,:])*kron((state - soln.hbar),(state - soln.hbar)))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state - soln.hbar),(state - soln.hbar)),(state - soln.hbar)))[1]
        end

        return y
    
    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R <: ThirdOrderSolutionStoch}

    function create_decision_rule(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)
        ny = length(soln.gbar)

        y = zeros(ny)
        
        for i = 1:ny
            y[i] = soln.gbar[i] + (soln.gx[i:i,:]*(state-soln.hbar))[1] + (1/2)*soln.gss[i] + (1/2)*((soln.gxx[i:i,:])*kron((state - soln.hbar),(state - soln.hbar)))[1] + (1/6)*soln.gsss[i]+ (3/6)*(soln.gssx[i:i,:]*(state - soln.hbar))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state - soln.hbar),(state - soln.hbar)),(state - soln.hbar)))[1]
        end

        return y
    
    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R <: Union{ChebyshevSolutionDet,ChebyshevSolutionStoch}}

    nx = length(soln.nodes)
    nv = length(soln.variables)
    ny = nv - nx

    T = eltype(soln.variables[1])

    w = Array{Array{T,nx},1}(undef,ny)
    if soln.node_generator == chebyshev_nodes
        for i = 1:ny
            w[i] = chebyshev_weights(soln.variables[nx+i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:ny
            w[i] = chebyshev_weights_extrema(soln.variables[nx+i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:ny
            w[i] = chebyshev_weights_extended(soln.variables[nx+i],soln.nodes,soln.order,soln.domain)
        end
    end

    function create_decision_rule(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end

        y = zeros(ny)
        
        for i = 1:ny
            y[i] = chebyshev_evaluate(w[i],state,soln.order,soln.domain)
        end

        return y
    
    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R <: Union{SmolyakSolutionDet,SmolyakSolutionStoch}}

    nx = size(soln.grid,2)
    nv = length(soln.variables)
    ny = nv - nx

    T = eltype(soln.variables[1])

    weights = Array{Array{T,1},1}(undef,ny)
    for i = 1:ny
        weights[i] = smolyak_weights(soln.variables[nx+i],soln.grid,soln.multi_index,soln.domain)
    end

    function create_decision_rule(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end

        y = zeros(ny)
        
        for i = 1:ny
            y[i] = smolyak_evaluate(weights[i],state,soln.multi_index,soln.domain)
        end

        return y
    
    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R <: Union{PiecewiseLinearSolutionDet,PiecewiseLinearSolutionStoch}}

    function create_decision_rule(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != length(soln.nodes)
            error("state vector has incorrect size")
        end
        nx = length(soln.nodes) 
        nv = length(soln.variables)
        ny = nv - nx

        y = zeros(ny)
        
        for i = 1:ny
            y[i] = piecewise_linear_evaluate(soln.variables[nx+i],soln.nodes,state)
        end

        return y
    
    end

    return create_decision_rule

end

function state_transition(soln::R) where {R <: FirstOrderSolutionDet}

    function create_state_transition(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        
        x_update = soln.hbar + soln.hx*(state-soln.hbar)

        return x_update
    
    end

    return create_state_transition

end

function state_transition(soln::R) where {R <: FirstOrderSolutionStoch}

    function create_state_transition(state::Array{T,1},shocks::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end

        x_update = soln.hbar + soln.hx*(state-soln.hbar) + soln.k*shocks

        return x_update
    
    end

    return create_state_transition

end

function state_transition(soln::R) where {R <: SecondOrderSolutionDet}

    function create_state_transition(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)

        x_update  = zeros(nx)
        
        for i = 1:nx
            x_update[i] = soln.hbar[i] + (soln.hx[i:i,:]*(state-soln.hbar))[1] + (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((state - soln.hbar),(state - soln.hbar)))
        end

        return x_update
    
    end

    return create_state_transition

end

function state_transition(soln::R) where {R <: SecondOrderSolutionStoch}

    function create_state_transition(state::Array{T,1},shocks::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end
        nx = length(soln.hbar)

        x_update  = zeros(nx)
        
        for i = 1:nx
            x_update[i] = soln.hbar[i] + (soln.hx[i:i,:]*(state-soln.hbar))[1] + (soln.k[i:i,:]*shocks)[1] + (1/2)*soln.hss[i] + (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((state - soln.hbar),(state - soln.hbar)))
        end

        return x_update
    
    end

    return create_state_transition

end

function state_transition(soln::R) where {R <: ThirdOrderSolutionDet}

    function create_state_transition(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)

        x_update = zeros(nx)
        
        for i = 1:nx
            x_update[i] = soln.hbar[i] + (soln.hx[i:i,:]*(state-soln.hbar))[1] + (1/2)*((soln.hxx[i:i,:])*kron((state - soln.hbar),(state - soln.hbar)))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state - soln.hbar),(state - soln.hbar)),(state - soln.hbar)))[1]
        end

        return x_update
    
    end

    return create_state_transition

end

function state_transition(soln::R) where {R <: ThirdOrderSolutionStoch}

    function create_state_transition(state::Array{T,1},shocks::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end
        nx = length(soln.hbar)

        x_update = zeros(nx)
        
        for i = 1:nx
            x_update[i] = soln.hbar[i] + (soln.hx[i:i,:]*(state-soln.hbar))[1] + (soln.k[i:i,:]*shocks)[1] + (1/2)*soln.hss[i] + (1/2)*((soln.hxx[i:i,:])*kron((state - soln.hbar),(state - soln.hbar)))[1] + (1/6)*soln.hsss[i] + (3/6)*(soln.hssx[i:i,:]*(state - soln.hbar))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state - soln.hbar),(state - soln.hbar)),(state - soln.hbar)))[1]
        end

        return x_update
    
    end

    return create_state_transition

end

function state_transition(soln::R) where {R <: ChebyshevSolutionDet}

    nx = length(soln.nodes)

    T = eltype(soln.variables[1])

    w = Array{Array{T,nx},1}(undef,nx)
    if soln.node_generator == chebyshev_nodes
        for i = 1:nx
            w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:nx
            w[i] = chebyshev_weights_extrema(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:nx
            w[i] = chebyshev_weights_extended(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    end

    function create_state_transition(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end

        x_update = zeros(nx)
        
        for i = 1:nx
            x_update[i] = chebyshev_evaluate(w[i],state,soln.order,soln.domain)
        end

        return x_update
    
    end

    return create_state_transition

end

function state_transition(soln::R) where {R <: ChebyshevSolutionStoch}

    nx = length(soln.nodes)
    ns = size(soln.k,2)

    T = eltype(soln.variables[1])

    w = Array{Array{T,nx},1}(undef,nx)
    if soln.node_generator == chebyshev_nodes
        for i = 1:nx
            w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:nx
            w[i] = chebyshev_weights_extrema(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:nx
            w[i] = chebyshev_weights_extended(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    end

    function create_state_transition(state::Array{T,1},shocks::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end

        x_update = zeros(nx)
        
        for i = 1:nx
            x_update[i] = chebyshev_evaluate(w[i],state,soln.order,soln.domain)
        end
        x_update[1:ns] += soln.k*shocks

        return x_update
    
    end

    return create_state_transition

end

function state_transition(soln::R) where {R <: SmolyakSolutionDet}

    nx = size(soln.grid,2)

    T = eltype(soln.variables[1])

    weights = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        weights[i] = smolyak_weights(soln.variables[i],soln.grid,soln.multi_index,soln.domain)
    end

    function create_state_transition(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end

        x_update = zeros(nx)
        
        for i = 1:nx
            x_update[i] = smolyak_evaluate(weights[i],state,soln.multi_index,soln.domain)
        end

        return x_update
    
    end

    return create_state_transition

end

function state_transition(soln::R) where {R <: SmolyakSolutionStoch}

    nx = size(soln.grid,2)
    ns = size(soln.k,2)

    T = eltype(soln.variables[1])

    weights = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        weights[i] = smolyak_weights(soln.variables[i],soln.grid,soln.multi_index,soln.domain)
    end

    function create_state_transition(state::Array{T,1},shocks::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end

        x_update = zeros(nx)
        
        for i = 1:nx
            x_update[i] = smolyak_evaluate(weights[i],state,soln.multi_index,soln.domain) 
        end
        x_update[1:ns] += soln.k*shocks

        return x_update
    
    end

    return create_state_transition

end

function state_transition(soln::R) where {R <: PiecewiseLinearSolutionDet}

    nx = length(soln.nodes) 

    function create_state_transition(state::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end

        x_update = zeros(nx)
        
        for i = 1:nx
            x_update[i] = piecewise_linear_evaluate(soln.variables[i],soln.nodes,state)
        end

        return x_update
    
    end

    return create_state_transition

end

function state_transition(soln::R) where {R <: PiecewiseLinearSolutionStoch}

    nx = length(soln.nodes) 
    ns = size(soln.k,2) 

    function create_state_transition(state::Array{T,1},shocks::Array{T,1}) where {T <: AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end

        x_update = zeros(nx)
        
        for i = 1:nx
            x_update[i] = piecewise_linear_evaluate(soln.variables[i],soln.nodes,state) 
        end
        x_update[1:ns] += soln.k*shocks

        return x_update
    
    end

    return create_state_transition

end

function expected_jumps(soln::R) where {R <: Union{PerturbationSolutionDet,ProjectionSolutionDet}}

    dec_rules   = decision_rule(soln)
    state_trans = state_transition(soln)

    function create_expected_jumps(state::Array{T,1}) where {T <: AbstractFloat}

        y = dec_rules(state_trans(state))

        return y
    
    end

    return create_expected_jumps

end

function expected_jumps(soln::R) where {R <: PerturbationSolutionStoch}

    dec_rules   = decision_rule(soln)
    state_trans = state_transition(soln)

    ns = size(soln.k,2) 

    function create_expected_jumps(state::Array{T,1}) where {T <: AbstractFloat}

        y = dec_rules(state_trans(state,zeros(ns)))

        return y
    
    end

    return create_expected_jumps

end

function expected_jumps(soln::R) where {R <: ChebyshevSolutionStoch}

    nx = length(soln.nodes)
    nv = length(soln.variables)
    ny = nv - nx
    ns = size(soln.k,2) 

    state_trans = state_transition(soln)

    T = eltype(soln.variables[1])

    weights        = Array{Array{T,nx},1}(undef,ny)
    scaled_weights = Array{Array{T,nx},1}(undef,ny)
    for i = 1:ny
        weights[i] = chebyshev_weights(soln.variables[nx+i],soln.nodes,soln.order,soln.domain)
    end

    if typeof(soln.integrals) == Array{T,ns}
        for i = 1:ny
            for j = 1:ns
                scaled_weights[i] = soln.integrals.*weights[i]
            end
        end
    else
        for i = 1:ny
            for j = 1:ns
                index = [1:ndims(weights[i]);]
                index[1],index[j] = index[j],index[1]
                scaled_weights[i] = permutedims(soln.integrals[j].*permutedims(weights[i],index),index)
            end
        end
    end

    function create_expected_jumps(state::Array{T,1}) where {T <: AbstractFloat}

        y = zeros(ny)

        for i = 1:ny

            y[i] = chebyshev_evaluate(scaled_weights[i],state_trans(state,zeros(ns)),soln.order,soln.domain)

        end

        return y

    end

    return create_expected_jumps

end

function expected_jumps(soln::R) where {R <: SmolyakSolutionStoch}

    nx = size(soln.grid,2)
    nv = length(soln.variables)
    ny = nv - nx
    ns = size(soln.k,2) 

    state_trans = state_transition(soln)

    T = eltype(soln.variables[1])

    weights = Array{Array{T,1},1}(undef,ny)
    for i = 1:ny
        weights[i] = smolyak_weights(soln.variables[nx+i],soln.grid,soln.multi_index,soln.domain).*soln.scale_factor
    end

    function create_expected_jumps(state::Array{T,1}) where {T <: AbstractFloat}

        y = zeros(ny)

        for i = 1:ny

            y[i] = smolyak_evaluate(weights[i],state_trans(state,zeros(ns)),soln.multi_index,soln.domain)

        end

        return y
    
    end

    return create_expected_jumps

end

function expected_jumps(soln::R) where {R <: PiecewiseLinearSolutionStoch}

    dec_rules   = decision_rule(soln)
    state_trans = state_transition(soln)

    ns = size(soln.k,2) 

    function create_expected_jumps(state::Array{T,1}) where {T <: AbstractFloat}

        y = dec_rules(state_trans(state,zeros(ns)))

        return y
    
    end

    return create_expected_jumps

end

function state_space_eqm(soln::R) where {R <: ModelSolution}

    dec_rules      = decision_rule(soln)
    transition_eqn = state_transition(soln)
    forecast_eqn   = expected_jumps(soln)

    eqm_dynamics = StateSpaceEqm(dec_rules,transition_eqn,forecast_eqn)

    return eqm_dynamics

end

function euler_errors(model::REModel,soln::R,domain::Union{Array{T,2},Array{T,1}},npoints::S,seed::S = 123456) where {S <: Integer, T <: AbstractFloat, R <: PerturbationSolutionDet}

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    dynamics = state_space_eqm(soln)

    euler_errors = zeros(length(model.eqns_approximated),npoints)

    Random.seed!(seed)
    states = domain[2,:] .+ rand(nx,npoints).*(domain[1,:] - domain[2,:])
    
    for i = 1:npoints
        state = states[:,i]
        point = [state;dynamics.g(state);dynamics.h(state);dynamics.gh(state)]
        f = model.dynamic_function(point)
        euler_errors[:,i] = f[model.eqns_approximated]
    end

    return euler_errors, states

end

function euler_errors(model::REModel,soln::R,domain::Union{Array{T,2},Array{T,1}},npoints::S,seed::S = 123456) where {S <: Integer, T <: AbstractFloat, R <: PerturbationSolutionStoch}

    nx = length(soln.hbar)
    ny = length(soln.gbar)
    ns = size(soln.k,2)

    shocks = zeros(ns)

    dynamics = state_space_eqm(soln)

    euler_errors = zeros(length(model.eqns_approximated),npoints)

    Random.seed!(seed)
    states = domain[2,:] .+ rand(nx,npoints).*(domain[1,:] - domain[2,:])
    
    for i = 1:npoints
        state = states[:,i]
        point = [state;dynamics.g(state);dynamics.h(state,shocks);dynamics.gh(state);shocks]
        f = model.dynamic_function(point)
        euler_errors[:,i] = f[model.eqns_approximated]
    end

    return euler_errors, states 

end

function euler_errors(model::REModel,soln::R,npoints::S,seed::S = 123456) where {S <: Integer, R <: ProjectionSolutionDet}

    nx = size(soln.domain,2)
    nv = length(soln.variables)
    ny = nv - nx

    dynamics = state_space_eqm(soln)

    euler_errors = zeros(length(model.eqns_approximated),npoints)

    Random.seed!(seed)
    states = soln.domain[2,:] .+ rand(nx,npoints).*(soln.domain[1,:] - soln.domain[2,:])
    
    for i = 1:npoints
        state = states[:,i]
        point = [state;dynamics.g(state);dynamics.h(state);dynamics.gh(state)]
        f = model.dynamic_function(point)
        euler_errors[:,i] = f[model.eqns_approximated]
    end

    return euler_errors, states

end

function euler_errors(model::REModel,soln::R,npoints::S,seed::S = 123456) where {S <: Integer, R <: ProjectionSolutionStoch}

    nx = size(soln.domain,2)
    nv = length(soln.variables)
    ny = nv - nx
    shocks = zeros(ns)

    dynamics = state_space_eqm(soln)

    euler_errors = zeros(length(model.eqns_approximated),npoints)

    Random.seed!(seed)
    states = soln.domain[2,:] .+ rand(nx,npoints).*(soln.domain[1,:] - soln.domain[2,:])
    
    for i = 1:npoints
        state = states[:,i]
        point = [state;dynamics.g(state);dynamics.h(state,shocks);dynamics.gh(state);shocks]
        f = model.dynamic_function(point)
        euler_errors[:,i] = f[model.eqns_approximated]
    end

    return euler_errors, states 

end
