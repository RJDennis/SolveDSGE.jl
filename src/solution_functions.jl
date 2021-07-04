# Solution functions

function compute_steady_state(model::REModel, x::Array{T,1}, tol::T, maxiters::S) where {T <: Real, S <: Integer}

    # Uses NLsolve (add this using the Package Manager if you don't have it) to
    # solve for the model's deterministic steady state

    if length(x) != model.number_equations
        error("The initialization has either too many or too few elements")
    end

    #equations = model.static_function
    #nlsoln = nlsolve(equations, x, xtol = tol, iterations = maxiters, autodiff = :forward, inplace = :false)
    equations = model.nlsolve_static_function
    nlsoln = nlsolve(equations, x, xtol = tol, iterations = maxiters, autodiff = :forward, inplace = :true)

    return nlsoln

end

##################################################################################

function compute_linearization(model::REModel, steady_state::Array{T,1}) where {T <: Real}

    # This function produces the Jacobian of the model's equations.

    equations = model.dynamic_function
    ns = model.number_shocks

    x = [steady_state;steady_state;zeros(ns)]
    d = ForwardDiff.jacobian(equations, x)

    return d

end

function solve_first_order(model::REModel,scheme::PerturbationScheme) # Follows Klein (2000)

    if scheme.order != "first"
        error("A first order perturbation must be specified")
    end

    ns = model.number_shocks
    if ns == 0
        soln = solve_first_order_det(model,scheme)
        return soln
    else
        soln = solve_first_order_stoch(model,scheme)
        return soln
    end

end

function solve_first_order_det(model::REModel,scheme::PerturbationScheme) # Follows Klein (2000)

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    steady_state = scheme.steady_state
    cutoff       = scheme.cutoff
    T            = eltype(steady_state)

    d = compute_linearization(model,steady_state)

    @views a =  d[1:nv,1:nv]
    @views b = -d[1:nv,nv+1:2*nv]

    r = schur(complex(a), complex(b))

    # Reorder the generalized eigenvalues so that those with modulus greater
    # than "cutoff" reside at the bottom.

    sel = (abs.(diag(r.S)./diag(r.T)).<cutoff)
    ordschur!(r, sel)

    s = r.S
    t = r.T
    #q = r.Q'  # So now q*a*z = s and q*b*z = t
    z = r.Z

    # Calculate the number of unstable eigenvalues

    grc = nv - sum(sel)

    soln_type = "determinate"
    if grc < ny
        soln_type = "indeterminate"
    elseif grc > ny
        soln_type = "explosive"
    end

    # Construct the rational expectations equilibrium by eliminating the
    # unstable dynamics

    @views z11 = z[1:nx,1:nx]
    @views z21 = z[nx+1:nv,1:nx]
    @views t11 = t[1:nx,1:nx]
    @views s11 = s[1:nx,1:nx]

    hx = real((z11/t11)*(s11/z11))
    gx = real(z21/z11)

    soln = FirstOrderSolutionDet(steady_state[1:nx],hx,steady_state[nx+1:nv],gx,grc,soln_type)
    return soln

end

function solve_first_order_stoch(model::REModel,scheme::PerturbationScheme) # Follows Klein (2000)

    nx = model.number_states
    ny = model.number_jumps
    ns = model.number_shocks
    nv = nx + ny

    steady_state = scheme.steady_state
    cutoff       = scheme.cutoff
    T            = eltype(steady_state)

    d = compute_linearization(model,steady_state)

    @views a =  d[1:nv,1:nv]
    @views b = -d[1:nv,nv+1:2*nv]

    r = schur(complex(a), complex(b))

    # Reorder the generalized eigenvalues so that those with modulus greater
    # than "cutoff" reside at the bottom.

    sel = (abs.(diag(r.S)./diag(r.T)).<cutoff)
    ordschur!(r, sel)

    s = r.S
    t = r.T
    #q = r.Q'  # So now q*a*z = s and q*b*z = t
    z = r.Z

    # Calculate the number of unstable eigenvalues

    grc = nv - sum(sel)

    soln_type = "determinate"
    if grc < ny
        soln_type = "indeterminate"
    elseif grc > ny
        soln_type = "explosive"
    end

    # Construct the rational expectations equilibrium by eliminating the
    # unstable dynamics

    @views z11 = z[1:nx,1:nx]
    @views z21 = z[nx+1:nv,1:nx]
    @views t11 = t[1:nx,1:nx]
    @views s11 = s[1:nx,1:nx]

    hx = real((z11/t11)*(s11/z11))
    gx = real(z21/z11)

    @views k = -d[1:nx,2*nv+1:2*nv+ns]
    sigma = eye(T,ns)
    soln = FirstOrderSolutionStoch(steady_state[1:nx],hx,k,steady_state[nx+1:nv],gx,sigma,grc,soln_type)
    return soln

end

function solve_second_order(model::REModel,scheme::PerturbationScheme) # Follows Gomme and Klein (2011)

    if scheme.order != "second"
        error("A second order perturbation must be specified")
    end

    ns = model.number_shocks
    if ns == 0
        soln = solve_second_order_det(model,scheme)
        return soln
    else
        soln = solve_second_order_stoch(model,scheme)
        return soln
    end

end

function solve_second_order_det(model::REModel,scheme::PerturbationScheme) # Follows Gomme and Klein (2011)

    nx = model.number_states
    ny = model.number_jumps
    ne = model.number_equations
    nv = nx + ny

    steady_state = scheme.steady_state
    cutoff       = scheme.cutoff
    T            = eltype(steady_state)

    first_order_soln = solve_first_order(model,PerturbationScheme(steady_state,cutoff,"first"))
    hx        = first_order_soln.hx
    gx        = first_order_soln.gx
    grc       = first_order_soln.grc
    soln_type = first_order_soln.soln_type

    # Calculate the first derivatives at the steady state

    d = compute_linearization(model,steady_state)

    # Calculate the Hessian at the steady state for each equation and stack the
    # Hessians vertically.

    point = [steady_state; steady_state]
    deriv2 = zeros(ne*2*nv,2*nv)
    for i = 1:nv
        deriv2[(i-1)*2*nv+1:i*2*nv,:] = ForwardDiff.hessian(model.each_eqn_function[i],point,ForwardDiff.HessianConfig(model.each_eqn_function[i],point,ForwardDiff.Chunk{1}()))[1:2*nv,1:2*nv]
     end

    # Construct partitioned first derivative matrices

    @views fx  = d[:,1:nx]
    @views fy  = d[:,(nx+1):nv]
    @views fxp = d[:,(nv+1):(nv+nx)]
    @views fyp = d[:,(nv+nx+1):2*nv]

    # Set up the Sylvester equation needed to construct the matrices on the
    # second-order states.

    m = [I; gx; hx; gx*hx]
    q = kron_prod_times_matrix(eye(T,nv),m',deriv2*m) #(eye(T,nv) ⊗ m')*(deriv2*m)
    b1 = kron(fxp,eye(T,nx)) #fxp ⊗ eye(T,nx) 
    b2 = kron(fyp,eye(T,nx)) #fyp ⊗ eye(T,nx)
    b4 = kron(fy,eye(T,nx)) #fy ⊗ eye(T,nx) 
    c1 = kron(eye(T,ny),hx') #eye(T,ny) ⊗ hx'
    c2 = kron(gx,eye(T,nx)) #gx ⊗ eye(T,nx)

    # Use a Sylvester equation solver to compute the second order terms on the
    # states.

    A = [b1+b2*c2 b4]
    C = hx
    B = [zeros(size(b2,1),nx*nx) b2*c1]
    D = q

    B .= A\B
    D .= A\D

    z = dsylvester(B,C,-D)
    hxx = z[1:nx*nx,:]
    gxx = z[nx^2+1:end,:]

    soln = SecondOrderSolutionDet(steady_state[1:nx],hx,hxx,steady_state[nx+1:nv],gx,gxx,grc,soln_type)
    return soln

end

function solve_second_order_stoch(model::REModel,scheme::PerturbationScheme) # Follows Gomme and Klein (2011)

    nx = model.number_states
    ny = model.number_jumps
    ns = model.number_shocks
    ne = model.number_equations
    nv = nx + ny

    steady_state = scheme.steady_state
    cutoff       = scheme.cutoff
    T            = eltype(steady_state)

    first_order_soln = solve_first_order(model,PerturbationScheme(steady_state,cutoff,"first"))
    hx        = first_order_soln.hx
    gx        = first_order_soln.gx
    grc       = first_order_soln.grc
    soln_type = first_order_soln.soln_type

    # Calculate the first derivatives at the steady state

    d = compute_linearization(model,steady_state)

    # Calculate the Hessian at the steady state for each equation and stack the
    # Hessians vertically.

    point = [steady_state; steady_state; zeros(ns)]
    deriv2 = zeros(ne*2*nv,2*nv)
    for i = 1:nv
        deriv2[(i-1)*2*nv+1:i*2*nv,:] = ForwardDiff.hessian(model.each_eqn_function[i],point,ForwardDiff.HessianConfig(model.each_eqn_function[i],point,ForwardDiff.Chunk{1}()))[1:2*nv,1:2*nv]
    end

    # Construct partitioned first derivative matrices

    @views fx  = d[:,1:nx]
    @views fy  = d[:,(nx+1):nv]
    @views fxp = d[:,(nv+1):(nv+nx)]
    @views fyp = d[:,(nv+nx+1):2*nv]

    # Set up the Sylvester equation needed to construct the matrices on the
    # second-order states.

    m = [I; gx; hx; gx*hx]
    q = kron_prod_times_matrix(eye(T,nv),m',deriv2*m) #(eye(T,nv) ⊗ m')*(deriv2*m)
    b1 = kron(fxp,eye(T,nx)) #fxp ⊗ eye(T,nx) 
    b2 = kron(fyp,eye(T,nx)) #fyp ⊗ eye(T,nx)
    b4 = kron(fy,eye(T,nx)) #fy ⊗ eye(T,nx) 
    c1 = kron(eye(T,ny),hx') #eye(T,ny) ⊗ hx'
    c2 = kron(gx,eye(T,nx)) #gx ⊗ eye(T,nx)

    # Use a Sylvester equation solver to compute the second order terms on the
    # states.

    A = [b1+b2*c2 b4]
    B = [zeros(size(b2,1),nx*nx) b2*c1]
    C = hx
    D = q

    B .= A\B
    D .= A\D

    z = dsylvester(B,C,-D)
    hxx = z[1:nx*nx,:]
    gxx = z[nx^2+1:end,:]

    # Set up the LP problem needed to construct the intercepts, which contain
    # the volatility effects

    k     = first_order_soln.k
    sigma = first_order_soln.sigma
    qq = [fxp+fyp*gx fyp+fy]
    nn = [zeros(nv,nx); I; gx]
    q = fyp*tracem(kron_prod_times_matrix(eye(T,ny),k*(sigma')*k',gxx)) + tracem(kron_prod_times_matrix(eye(T,nv),nn',deriv2*nn*k*(sigma')*k'))
    ss = -qq\q
    hss = ss[1:nx]
    gss = ss[nx+1:nv]

    soln = SecondOrderSolutionStoch(steady_state[1:nx],hx,hss,hxx,k,steady_state[nx+1:nv],gx,gss,gxx,sigma,grc,soln_type)
    return soln

end

function solve_third_order(model::REModel,scheme::PerturbationScheme, skewness::Union{Array{T,1},Array{T,2}} = zeros(model.number_shocks,model.number_shocks^2)) where {T <: Real} # Follows Binning (2013)

    if scheme.order != "third"
        error("A third order perturbation must be supplied")
    end

    ns = model.number_shocks
    if ns == 0
        soln = solve_third_order_det(model,scheme)
        return soln
    else
        soln = solve_third_order_stoch(model,scheme,skewness)
        return soln
    end

end

function solve_third_order_det(model::REModel, scheme::PerturbationScheme) # Follows Binning (2013)

    nv = model.number_variables
    nx = model.number_states
    ny = model.number_jumps
    ne = model.number_equations

    steady_state = scheme.steady_state
    cutoff       = scheme.cutoff

    first_order_soln = solve_first_order(model,PerturbationScheme(steady_state,cutoff,"first"))
    hx        = first_order_soln.hx
    gx        = first_order_soln.gx
    grc       = first_order_soln.grc
    soln_type = first_order_soln.soln_type

    # Compute the first, second, and third derivatives

    model_equations = model.each_eqn_function
    point = [steady_state; steady_state]

    first_derivs  = zeros(ne,2*nv)
    second_derivs = zeros(ne,4*nv^2)
    third_derivs  = zeros(ne,8*nv^3)

    for i = 1:ne

        first_d(x) = ForwardDiff.gradient(model_equations[i],x,ForwardDiff.GradientConfig(model_equations[i],x,ForwardDiff.Chunk{1}()))[1:2*nv]
        first_derivs[i,:] = Matrix(first_d(point)')

        second_d(x) = ForwardDiff.hessian(model_equations[i],x,ForwardDiff.HessianConfig(model_equations[i],x,ForwardDiff.Chunk{1}()))[1:2*nv,1:2*nv]
        #second_d(x) = ForwardDiff.jacobian(first_d,x,ForwardDiff.JacobianConfig(first_d,x,ForwardDiff.Chunk{1}()))[1:2*n,1:2*n]
        second_derivs[i,:] = Matrix(vec(second_d(point))')

        third_d(x) = ForwardDiff.jacobian(second_d,x,ForwardDiff.JacobianConfig(second_d,x,ForwardDiff.Chunk{1}()))[:,1:2*nv]
        third_derivs[i,:] = Matrix(vec(third_d(point))')

    end

    # Solve for the terms in the second-order solution

    # Construct the Mx matrix --- Gradient of the policy functions

    Mx = [I; gx; hx; gx*hx]

    @views fx  = first_derivs[:,1:nx]
    @views fy  = first_derivs[:,nx+1:nv]
    @views fxp = first_derivs[:,nv+1:nv+nx]
    @views fyp = first_derivs[:,nv+nx+1:2*nv]

    # Compute the second-order terms on the state variables

    A = [fxp+fyp*gx fy]
    B = [zeros(nv,nx) fyp]
    D = -matrix_times_kron_prod(second_derivs,Mx,Mx)

    z = martin_van_loan(A,B,hx,D,1)
    hxx = z[1:nx,:]
    gxx = z[nx+1:nv,:]

    # Construct the Mxx matrix --- Hessian of the policy functions

    Mxx = [zeros(nx,nx^2); gxx; hxx; gxx*kron(hx,hx) + gx*hxx]

    # Construct the hxx_dagger matrix

    hxx_dagger = zeros(nx^2,nx^3)
    for i = 1:nx
        @views hxx_dagger[:,(i-1)*nx^2+1:i*nx^2] = kron(hx,hxx[:,(i-1)*nx+1:i*nx])
    end

    # Solve for the third order coefficients on the state variables

    omega = create_omega3(nx)

    #A = [fxp+fyp*gx fy]
    #B = [zeros(nv,nx) fyp]
    D = (-matrix_times_kron_prod(third_derivs,Mx,Mx,Mx) - matrix_times_kron_prod(second_derivs,Mx,Mxx)*omega - fyp*gxx*(kron(hx,hxx) + kron(hxx,hx) + hxx_dagger))

    z = martin_van_loan(A,B,hx,D,2)
    hxxx = z[1:nx,:]
    gxxx = z[nx+1:nv,:]

    soln = ThirdOrderSolutionDet(steady_state[1:nx],hx,hxx,hxxx,steady_state[nx+1:nv],gx,gxx,gxxx,grc,soln_type)
    return soln

end

function solve_third_order_stoch(model::REModel, scheme::PerturbationScheme, skewness::Union{Array{T,1},Array{T,2}}) where {T <: Real} # Follows Binning (2013)

    ns = model.number_shocks
    nv = model.number_variables
    nx = model.number_states
    ny = model.number_jumps
    ne = model.number_equations

    steady_state = scheme.steady_state
    cutoff       = scheme.cutoff

    first_order_soln = solve_first_order(model,PerturbationScheme(steady_state,cutoff,"first"))
    hx        = first_order_soln.hx
    gx        = first_order_soln.gx
    grc       = first_order_soln.grc
    soln_type = first_order_soln.soln_type

    # Compute the first, second, and third derivatives

    model_equations = model.each_eqn_function
    point = [steady_state; steady_state; zeros(ns)]

    first_derivs  = zeros(ne,2*nv)
    second_derivs = zeros(ne,4*nv^2)
    third_derivs  = zeros(ne,8*nv^3)

    for i = 1:ne

        first_d(x) = ForwardDiff.gradient(model_equations[i],x,ForwardDiff.GradientConfig(model_equations[i],x,ForwardDiff.Chunk{1}()))[1:2*nv]
        first_derivs[i,:] = Matrix(first_d(point)')

        second_d(x) = ForwardDiff.hessian(model_equations[i],x,ForwardDiff.HessianConfig(model_equations[i],x,ForwardDiff.Chunk{1}()))[1:2*nv,1:2*nv]
        #second_d(x) = ForwardDiff.jacobian(first_d,x,ForwardDiff.JacobianConfig(first_d,x,ForwardDiff.Chunk{1}()))[1:2*n,1:2*n]
        second_derivs[i,:] = Matrix(vec(second_d(point))')

        third_d(x) = ForwardDiff.jacobian(second_d,x,ForwardDiff.JacobianConfig(second_d,x,ForwardDiff.Chunk{1}()))[:,1:2*nv]
        third_derivs[i,:] = Matrix(vec(third_d(point))')

    end

    # Solve for the terms in the second-order solution

    # Construct the Mx matrix --- Gradient of the policy functions

    Mx = [I; gx; hx; gx*hx]

    @views fx  = first_derivs[:,1:nx]
    @views fy  = first_derivs[:,nx+1:nv]
    @views fxp = first_derivs[:,nv+1:nv+nx]
    @views fyp = first_derivs[:,nv+nx+1:2*nv]

    # Compute the second-order terms on the state variables

    A = [fxp+fyp*gx fy]
    B = [zeros(nv,nx) fyp]
    D = -matrix_times_kron_prod(second_derivs,Mx,Mx)

    z = martin_van_loan(A,B,hx,D,1)
    hxx = z[1:nx,:]
    gxx = z[nx+1:nv,:]

    # Construct the Mxx matrix --- Hessian of the policy functions

    Mxx = [zeros(nx,nx^2); gxx; hxx; gxx*kron(hx,hx) + gx*hxx]

    # Construct the hxx_dagger matrix

    hxx_dagger = zeros(nx^2,nx^3)
    for i = 1:nx
        @views hxx_dagger[:,(i-1)*nx^2+1:i*nx^2] .= kron(hx,hxx[:,(i-1)*nx+1:i*nx])
    end

    # Solve for the third order coefficients on the state variables

    omega = create_omega3(nx)

    #A = [fxp+fyp*gx fy]
    #B = [zeros(nv,nx) fyp]
    D = (-matrix_times_kron_prod(third_derivs,Mx,Mx,Mx) .- matrix_times_kron_prod(second_derivs,Mx,Mxx)*omega .- fyp*gxx*(kron(hx,hxx) .+ kron(hxx,hx) .+ hxx_dagger))

    z = martin_van_loan(A,B,hx,D,2)
    hxxx = z[1:nx,:]
    gxxx = z[nx+1:nv,:]

    # Compute the risk-adjustment terms

    k     = first_order_soln.k
    sigma = first_order_soln.sigma

    Ns = [zeros(nv,nx); I; gx]

    A = [fxp+fyp*gx fy+fyp]
    F = -trm(matrix_times_kron_prod(second_derivs,Ns,Ns*k*sigma*k')) .- fyp*trm(matrix_times_kron_prod(gxx,eye(T,nx),k*sigma*k'))

    z = A\F
    hss = reshape(z[1:nx,:],nx)
    gss = reshape(z[nx+1:nv,:],ny)

    # Solve for the volatility-related terms.

    Nsx = [zeros(nv+nx,nx^2); gxx*kron(hx,eye(T,nx))]
    Pss = [zeros(nx,1); gss; hss; gx*hss + trm(gxx*kron(eye(T,nx),k*sigma*k')) + gss]

    C = hx
    F = (-trm2(matrix_times_kron_prod(third_derivs,Mx,Ns,Ns*k*sigma*k')) - 2*trm2(matrix_times_kron_prod(second_derivs,Nsx,Ns*k*sigma*k')) - matrix_times_kron_prod(second_derivs,Mx,Pss[:,:])
      -fyp*(trm2(gxxx*kron(hx,eye(T,nx^2))*kron(eye(T,nx^2),k*sigma*k')) + matrix_times_kron_prod(gxx,hx,hss[:,:])))

    F = A\F

    z = dsylvester(B,C,F)
    hssx = z[1:nx,:]
    gssx = z[nx+1:nv,:]

    # Now check to see whether skewness is relevant.

    hsss = zeros(nx)
    gsss = zeros(ny)

    if sum(isequal.(skewness,0.0)) == 0

        skew = zeros(nx,nx^2)
        skew[1:ns,1:ns^2] = skewness

        Nss = [zeros(nv+nx,nx^2); gxx]

        #A = [fxp+fyp*gx fy+fyp]
        F = -trm(third_derivs*kron(Ns,kron(Ns,Ns*skew))) - 3*trm(second_derivs*kron(Nss,Ns*skew)) - fyp*trm(gxxx*kron(eye(T,nx^2),skew))
        z = A\F
        hsss = reshape(z[1:nx,:],nx)
        gsss = reshape(z[nx+1:nv,:],ny)

    end

    soln = ThirdOrderSolutionStoch(steady_state[1:nx],hx,hss,hxx,hsss,hssx,hxxx,k,steady_state[nx+1:nv],gx,gss,gxx,gsss,gssx,gxxx,sigma,skewness,grc,soln_type)
    return soln

end

function solve_nonlinear(model::REModel,scheme::ChebyshevSchemeDet)

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number    = scheme.node_number
    order          = scheme.order
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    state = Array{T,1}(undef,nx)

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],Tuple(length.(grid)))
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],Tuple(length.(grid)))
    end

    weights = Array{Array{T,nx},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    init  = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        if node_generator == chebyshev_nodes
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extrema
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extrema(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extended
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extended(variables[jumps_approximated[i]],grid,order,domain)
            end
        end

        for i = 1:N

            ind = ind2sub(i,Tuple(length.(grid)))
            for j = 1:nx
                state[j] = grid[j][ind[j]]
            end
            for j = 1:nv
                init[j] = variables[j][ind...]
            end

            projection_equations = model.closure_function(state,weights,order,domain,chebyshev_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][ind...] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionDet([variables[ny+1:end];variables[1:ny]],weights,grid,order,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::ChebyshevSchemeDet,threads::S) where {S <: Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number    = scheme.node_number
    order          = scheme.order
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],Tuple(length.(grid)))
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],Tuple(length.(grid)))
    end

    weights = Array{Array{T,nx},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        if node_generator == chebyshev_nodes
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extrema
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extrema(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extended
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extended(variables[jumps_approximated[i]],grid,order,domain)
            end
        end

        @sync @qthreads for t = 1:threads
            for i = t:threads:N

                sub = ind2sub(i,Tuple(length.(grid)))
                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init  = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function(state,weights,order,domain,chebyshev_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionDet([variables[ny+1:end];variables[1:ny]],weights,grid,order,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::ChebyshevSchemeStoch)

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number    = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    order          = scheme.order
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    d = compute_linearization(model,initial_guess)
    k = -d[1:ns,2*nv+1:end]
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = Array{Array{T,1},1}(undef,ns)
    for i = 1:ns
        if length(order) == 1
            integrals[i] = compute_chebyshev_integrals(eps_nodes,eps_weights,grid[i],order,-d[i,i],k[i,i])
        else
            integrals[i] = compute_chebyshev_integrals(eps_nodes,eps_weights,grid[i],order[i],-d[i,i],k[i,i])
        end
    end

    state = Array{T,1}(undef,nx)

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],Tuple(length.(grid)))
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],Tuple(length.(grid)))
    end

    weights        = Array{Array{T,nx},1}(undef,length(jumps_approximated))
    scaled_weights = Array{Array{T,nx},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    init  = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        if node_generator == chebyshev_nodes
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extrema
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extrema(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extended
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extended(variables[jumps_approximated[i]],grid,order,domain)
            end
        end

        for i = 1:length(jumps_approximated)
            for j = 1:ns
                index = [1:ndims(weights[i]);]
                index[1],index[j] = index[j],index[1]
                scaled_weights[i] = permutedims(integrals[j].*permutedims(weights[i],index),index)
            end
        end

        for i = 1:N

            ind = ind2sub(i,Tuple(length.(grid)))
            for j = 1:nx
                state[j] = grid[j][ind[j]]
            end
            for j = 1:nv
                init[j] = variables[j][ind...]
            end

            projection_equations = model.closure_function(state,scaled_weights,order,domain,chebyshev_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionStoch([variables[ny+1:end];variables[1:ny]],weights,grid,order,domain,Matrix(k*k'),iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::ChebyshevSchemeStoch,threads::S) where {S <: Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number    = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    order          = scheme.order
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    d = compute_linearization(model,initial_guess)
    k = -d[1:ns,2*nv+1:end]
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = Array{Array{T,1},1}(undef,ns)
    for i = 1:ns
        if length(order) == 1
            integrals[i] = compute_chebyshev_integrals(eps_nodes,eps_weights,grid[i],order,-d[i,i],k[i,i])
        else
            integrals[i] = compute_chebyshev_integrals(eps_nodes,eps_weights,grid[i],order[i],-d[i,i],k[i,i])
        end
    end

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],Tuple(length.(grid)))
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],Tuple(length.(grid)))
    end

    weights        = Array{Array{T,nx},1}(undef,length(jumps_approximated))
    scaled_weights = Array{Array{T,nx},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        if node_generator == chebyshev_nodes
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extrema
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extrema(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extended
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extended(variables[jumps_approximated[i]],grid,order,domain)
            end
        end

        for i = 1:length(jumps_approximated)
            for j = 1:ns
                index = [1:ndims(weights[i]);]
                index[1],index[j] = index[j],index[1]
                scaled_weights[i] = permutedims(integrals[j].*permutedims(weights[i],index),index)
            end
        end

        @sync @qthreads for t = 1:threads
            for i = t:threads:N

                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end

                init  = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function(state,scaled_weights,order,domain,chebyshev_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionStoch([variables[ny+1:end];variables[1:ny]],weights,grid,order,domain,Matrix(k*k'),iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::ChebyshevSchemeDet) where {R <: PerturbationSolutionDet}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    node_generator = scheme.node_generator
    node_number    = scheme.node_number
    order          = scheme.order
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    hbar = soln.hbar
    hx   = soln.hx
    gbar = soln.gbar
    gx   = soln.gx

    T = typeof(scheme.tol_fix_point_solver)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    state = Array{T,1}(undef,nx)

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    for j = 1:N
        sub = ind2sub(j,Tuple(length.(grid)))

        for i = 1:nx
            state[i] = grid[i][sub[i]]
        end

        for i = 1:ny
            variables[i][j] = gbar[i] + (gx[i:i,:]*(state - hbar))[1]
        end
        for i = 1:nx
            variables[ny+i][j] = hbar[i] + (hx[i:i,:]*(state - hbar))[1]
        end

        if typeof(soln) <: SecondOrderSolutionDet
            for i = 1:ny
                variables[i][j] += (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
            end
        end

        if typeof(soln) <: ThirdOrderSolutionDet
            for i = 1:ny
                variables[i][j] += (1/2)*((soln.gxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*((soln.hxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
            end
        end

    end

    weights = Array{Array{T,nx},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    init  = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        if node_generator == chebyshev_nodes
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extrema
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extrema(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extended
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extended(variables[jumps_approximated[i]],grid,order,domain)
            end
        end

        for i = 1:N

            sub = ind2sub(i,Tuple(length.(grid)))
            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function(state,weights,order,domain,chebyshev_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionDet([variables[ny+1:end];variables[1:ny]],weights,grid,order,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::ChebyshevSchemeDet,threads::S) where {R <: PerturbationSolutionDet, S <: Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    node_generator = scheme.node_generator
    node_number    = scheme.node_number
    order          = scheme.order
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    hbar = soln.hbar
    hx   = soln.hx
    gbar = soln.gbar
    gx   = soln.gx

    T = typeof(scheme.tol_fix_point_solver)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    @sync @qthreads for t = 1:threads
        for j = t:threads:N
            sub = ind2sub(j,Tuple(length.(grid)))

            state = Array{T,1}(undef,nx)
            for i = 1:nx
                state[i] = grid[i][sub[i]]
            end

            for i = 1:ny
                variables[i][j] = gbar[i] + (gx[i:i,:]*(state - hbar))[1]
            end
            for i = 1:nx
                variables[ny+i][j] = hbar[i] + (hx[i:i,:]*(state - hbar))[1]
            end

            if typeof(soln) <: SecondOrderSolutionDet
                for i = 1:ny
                    variables[i][j] += (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
                end
                for i = 1:nx
                    variables[ny+i][j] += (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
                end
            end

            if typeof(soln) <: ThirdOrderSolutionDet
                for i = 1:ny
                    variables[i][j] += (1/2)*((soln.gxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
                end
                for i = 1:nx
                    variables[ny+i][j] += (1/2)*((soln.hxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
                end
            end
        end
    end

    weights = Array{Array{T,nx},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        if node_generator == chebyshev_nodes
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extrema
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extrema(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extended
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extended(variables[jumps_approximated[i]],grid,order,domain)
            end
        end

        @sync @qthreads for t = 1:threads
            for i = t:threads:N

                sub = ind2sub(i,Tuple(length.(grid)))
                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init  = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function(state,weights,order,domain,chebyshev_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionDet([variables[ny+1:end];variables[1:ny]],weights,grid,order,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::ChebyshevSchemeStoch) where {R <: PerturbationSolutionStoch}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    node_generator = scheme.node_generator
    node_number    = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    order          = scheme.order
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    hbar = soln.hbar
    hx   = soln.hx
    k    = soln.k
    gbar = soln.gbar
    gx   = soln.gx

    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        if typeof(soln) <: FirstOrderSolutionStoch
            soln_fo = soln
        else
            soln_fo = FirstOrderSolutionStoch(hbar,hx,k,gbar,gx,soln.sigma,soln.grc,soln.soln_type)
        end

        state_vars, jump_vars = compute_variances(soln_fo)

        domain = Matrix([hbar + 3*sqrt.(diag(state_vars)) hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = Array{Array{T,1},1}(undef,ns)
    for i = 1:ns
        if length(order) == 1
            integrals[i] = compute_chebyshev_integrals(eps_nodes,eps_weights,grid[i],order,hx[i,i],k[i,i])
        else
            integrals[i] = compute_chebyshev_integrals(eps_nodes,eps_weights,grid[i],order[i],hx[i,i],k[i,i])
        end
    end

    state = Array{T,1}(undef,nx)

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    for j = 1:N
        sub = ind2sub(j,Tuple(length.(grid)))

        for i = 1:nx
            state[i] = grid[i][sub[i]]
        end

        for i = 1:ny
            variables[i][j] = gbar[i] + (gx[i:i,:]*(state - hbar))[1]
        end
        for i = 1:nx
            variables[ny+i][j] = hbar[i] + (hx[i:i,:]*(state - hbar))[1]
        end

        if typeof(soln) <: SecondOrderSolutionStoch
            for i = 1:ny
                variables[i][j] += (1/2)*soln.gss[i] + (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*soln.hss[i] + (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
            end
        end

        if typeof(soln) <: ThirdOrderSolutionStoch
            for i = 1:ny
                variables[i][j] += (1/2)*soln.gss[i] + (1/2)*((soln.gxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (3/6)*(soln.gssx[i:i,:]*(state - hbar))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*soln.hss[i] + (1/2)*((soln.hxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (3/6)*(soln.hssx[i:i,:]*(state - hbar))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
            end
        end

    end

    weights        = Array{Array{T,nx},1}(undef,length(jumps_approximated))
    scaled_weights = Array{Array{T,nx},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    init  = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        if node_generator == chebyshev_nodes
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extrema
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extrema(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extended
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extended(variables[jumps_approximated[i]],grid,order,domain)
            end
        end

        for i = 1:length(jumps_approximated)
            for j = 1:ns
                index = [1:ndims(weights[i]);]
                index[1],index[j] = index[j],index[1]
                scaled_weights[i] = permutedims(integrals[j].*permutedims(weights[i],index),index)
            end
        end

        for i = 1:N

            sub = ind2sub(i,Tuple(length.(grid)))
            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function(state,scaled_weights,order,domain,chebyshev_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionStoch([variables[ny+1:end];variables[1:ny]],weights,grid,order,domain,Matrix(k[1:ns,:]*k[1:ns,:]'),iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::ChebyshevSchemeStoch,threads::S) where {R <: PerturbationSolutionStoch, S <: Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    node_generator = scheme.node_generator
    node_number    = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    order          = scheme.order
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    hbar = soln.hbar
    hx   = soln.hx
    k    = soln.k
    gbar = soln.gbar
    gx   = soln.gx

    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        if typeof(soln) <: FirstOrderSolutionStoch
            soln_fo = soln
        else
            soln_fo = FirstOrderSolutionStoch(hbar,hx,k,gbar,gx,soln.sigma,soln.grc,soln.soln_type)
        end

        state_vars, jump_vars = compute_variances(soln_fo)

        domain = Matrix([hbar + 3*sqrt.(diag(state_vars)) hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = Array{Array{T,1},1}(undef,ns)
    for i = 1:ns
        if length(order) == 1
            integrals[i] = compute_chebyshev_integrals(eps_nodes,eps_weights,grid[i],order,hx[i,i],k[i,i])
        else
            integrals[i] = compute_chebyshev_integrals(eps_nodes,eps_weights,grid[i],order[i],hx[i,i],k[i,i])
        end
    end

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    for j = 1:N

        sub = ind2sub(j,Tuple(length.(grid)))

        state = Array{T,1}(undef,nx)
        for i = 1:nx
            state[i] = grid[i][sub[i]]
        end

        for i = 1:ny
            variables[i][j] = gbar[i] + (gx[i:i,:]*(state - hbar))[1]
        end
        for i = 1:nx
            variables[ny+i][j] = hbar[i] + (hx[i:i,:]*(state - hbar))[1]
        end

        if typeof(soln) <: SecondOrderSolutionStoch
            for i = 1:ny
                variables[i][j] += (1/2)*soln.gss[i] + (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*soln.hss[i] + (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
            end
        end

        if typeof(soln) <: ThirdOrderSolutionStoch
            for i = 1:ny
                variables[i][j] += (1/2)*soln.gss[i] + (1/2)*((soln.gxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (3/6)*(soln.gssx[i:i,:]*(state - hbar))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*soln.hss[i] + (1/2)*((soln.hxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (3/6)*(soln.hssx[i:i,:]*(state - hbar))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
            end
        end
    end

    weights        = Array{Array{T,nx},1}(undef,length(jumps_approximated))
    scaled_weights = Array{Array{T,nx},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        if node_generator == chebyshev_nodes
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extrema
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extrema(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extended
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extended(variables[jumps_approximated[i]],grid,order,domain)
            end
        end

        for i = 1:length(jumps_approximated)
            for j = 1:ns
                index = [1:ndims(weights[i]);]
                index[1],index[j] = index[j],index[1]
                scaled_weights[i] = permutedims(integrals[j].*permutedims(weights[i],index),index)
            end
        end

        @sync @qthreads for t = 1:threads
            for i = t:threads:N

                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end

                init  = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function(state,scaled_weights,order,domain,chebyshev_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionStoch([variables[ny+1:end];variables[1:ny]],weights,grid,order,domain,Matrix(k[1:ns,:]*k[1:ns,:]'),iters,scheme.node_generator)


    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::ChebyshevSchemeDet) where {R <: ProjectionSolutionDet}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number    = scheme.node_number
    order          = scheme.order
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        domain = soln.domain
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    state = Array{T,1}(undef,nx)

    N = prod(length.(grid))

    soln_variables = [soln.variables[nx+1:end];soln.variables[1:nx]]

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
        if typeof(soln) <: ChebyshevSolutionDet
            if soln.node_generator == chebyshev_nodes
                w = chebyshev_weights(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extrema
                w = chebyshev_weights_extrema(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extended
                w = chebyshev_weights_extended(soln_variables[i],soln.nodes,soln.order,soln.domain)
            end

            for j = 1:N
                sub = ind2sub(j,Tuple(length.(grid)))

                for l = 1:nx
                    state[l] = grid[l][sub[l]]
                end

                variables[i][j] = chebyshev_evaluate(w,state,soln.order,domain)

            end
        elseif typeof(soln) <: SmolyakSolutionDet
            w = smolyak_weights(soln_variables[i],soln.grid,soln.multi_index,soln.domain)

            for j = 1:N
                sub = ind2sub(j,Tuple(length.(grid)))

                for l = 1:nx
                    state[l] = grid[l][sub[l]]
                end

                variables[i][j] = smolyak_evaluate(w,state,soln.multi_index,domain)

            end

        elseif typeof(soln) <: PiecewiseLinearSolutionDet

            for j = 1:N
                sub = ind2sub(j,Tuple(length.(grid)))

                for l = 1:nx
                    state[l] = grid[l][sub[l]]
                end

                variables[i][j] = piecewise_linear_evaluate(soln_variables[i],soln.nodes,state)

            end

        end

    end

    weights = Array{Array{T,nx},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    init  = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        if node_generator == chebyshev_nodes
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extrema
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extrema(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extended
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extended(variables[jumps_approximated[i]],grid,order,domain)
            end
        end

        for i = 1:N

            sub = ind2sub(i,Tuple(length.(grid)))
            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function(state,weights,order,domain,chebyshev_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end

        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionDet([variables[ny+1:end];variables[1:ny]],weights,grid,order,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::ChebyshevSchemeDet,threads::S) where {R <: ProjectionSolutionDet, S <: Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number    = scheme.node_number
    order          = scheme.order
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        domain = soln.domain
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    N = prod(length.(grid))

    soln_variables = [soln.variables[nx+1:end];soln.variables[1:nx]]

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
        if typeof(soln) <: ChebyshevSolutionDet
            if soln.node_generator == chebyshev_nodes
                w = chebyshev_weights_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extrema
                w = chebyshev_weights_extrema_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extended
                w = chebyshev_weights_extended_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            end

            @sync @qthreads for t = 1:threads
                for j = t:threads:N
                    sub = ind2sub(j,Tuple(length.(grid)))

                    state = Array{T,1}(undef,nx)
                    for l = 1:nx
                        state[l] = grid[l][sub[l]]
                    end

                    variables[i][j] = chebyshev_evaluate(w,state,soln.order,domain)

                end
            end
        elseif typeof(soln) <: SmolyakSolutionDet
            w = smolyak_weights_threaded(soln_variables[i],soln.grid,soln.multi_index,soln.domain)

            @sync @qthreads for t = 1:threads
                for j = t:threads:N
                    sub = ind2sub(j,Tuple(length.(grid)))

                    state = Array{T,1}(undef,nx)
                    for l = 1:nx
                        state[l] = grid[l][sub[l]]
                    end

                    variables[i][j] = smolyak_evaluate(w,state,soln.multi_index,domain)

                end
            end

        elseif typeof(soln) <: PiecewiseLinearSolutionDet

            @sync @qthreads for t = 1:threads
                for j = t:threads:N
                    sub = ind2sub(j,Tuple(length.(grid)))

                    state = Array{T,1}(undef,nx)
                    for l = 1:nx
                        state[l] = grid[l][sub[l]]
                    end

                    variables[i][j] = piecewise_linear_evaluate(soln_variables[i],soln.nodes,state)

                end
            end
        end
    end

    weights = Array{Array{T,nx},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        if node_generator == chebyshev_nodes
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extrema
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extrema(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extended
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extended(variables[jumps_approximated[i]],grid,order,domain)
            end
        end

        @sync @qthreads for t = 1:threads
            for i = t:threads:N

                sub = ind2sub(i,Tuple(length.(grid)))
                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init  = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function(state,weights,order,domain,chebyshev_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end

            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionDet([variables[ny+1:end];variables[1:ny]],weights,grid,order,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::ChebyshevSchemeStoch) where {R <: ProjectionSolutionStoch}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number    = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    order          = scheme.order
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        domain = soln.domain
    end

    d = compute_linearization(model,initial_guess)
    k = -d[1:ns,2*nv+1:end]
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = Array{Array{T,1},1}(undef,ns)
    for i = 1:ns
        if length(order) == 1
            integrals[i] = compute_chebyshev_integrals(eps_nodes,eps_weights,grid[i],order,-d[i,i],k[i,i])
        else
            integrals[i] = compute_chebyshev_integrals(eps_nodes,eps_weights,grid[i],order[i],-d[i,i],k[i,i])
        end
    end

    state = Array{T,1}(undef,nx)

    N = prod(length.(grid))

    soln_variables = [soln.variables[nx+1:end];soln.variables[1:nx]]

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
        if typeof(soln) <: ChebyshevSolutionStoch
            if soln.node_generator == chebyshev_nodes
                w = chebyshev_weights(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extrema
                w = chebyshev_weights_extrema(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extended
                w = chebyshev_weights_extended(soln_variables[i],soln.nodes,soln.order,soln.domain)
            end

            for j = 1:N
                sub = ind2sub(j,Tuple(length.(grid)))

                for l = 1:nx
                    state[l] = grid[l][sub[l]]
                end

                variables[i][j] = chebyshev_evaluate(w,state,soln.order,domain)

            end
        elseif typeof(soln) <: SmolyakSolutionStoch
            w = smolyak_weights(soln_variables[i],soln.grid,soln.multi_index,soln.domain)

            for j = 1:N
                sub = ind2sub(j,Tuple(length.(grid)))

                for l = 1:nx
                    state[l] = grid[l][sub[l]]
                end

                variables[i][j] = smolyak_evaluate(w,state,soln.multi_index,domain)

            end

        elseif typeof(soln) <: PiecewiseLinearSolutionStoch

            for j = 1:N
                sub = ind2sub(j,Tuple(length.(grid)))

                for l = 1:nx
                    state[l] = grid[l][sub[l]]
                end

                variables[i][j] = piecewise_linear_evaluate(soln_variables[i],soln.nodes,state)

            end

        end

    end

    weights        = Array{Array{T,nx},1}(undef,length(jumps_approximated))
    scaled_weights = Array{Array{T,nx},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    init  = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        if node_generator == chebyshev_nodes
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extrema
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extrema(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extended
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extended(variables[jumps_approximated[i]],grid,order,domain)
            end
        end

        for i = 1:length(jumps_approximated)
            for j = 1:ns
                index = [1:ndims(weights[i]);]
                index[1],index[j] = index[j],index[1]
                scaled_weights[i] = permutedims(integrals[j].*permutedims(weights[i],index),index)
            end
        end

        for i = 1:N
            sub = ind2sub(i,Tuple(length.(grid)))
            
            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function(state,scaled_weights,order,domain,chebyshev_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end

        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionStoch([variables[ny+1:end];variables[1:ny]],weights,grid,order,domain,Matrix(k*k'),iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::ChebyshevSchemeStoch,threads::S) where {R <: ProjectionSolutionStoch, S <: Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number    = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    order          = scheme.order
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        domain = soln.domain
    end

    d = compute_linearization(model,initial_guess)
    k = -d[1:ns,2*nv+1:end]
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = Array{Array{T,1},1}(undef,ns)
    for i = 1:ns
        if length(order) == 1
            integrals[i] = compute_chebyshev_integrals(eps_nodes,eps_weights,grid[i],order,-d[i,i],k[i,i])
        else
            integrals[i] = compute_chebyshev_integrals(eps_nodes,eps_weights,grid[i],order[i],-d[i,i],k[i,i])
        end
    end

    N = prod(length.(grid))

    soln_variables = [soln.variables[nx+1:end];soln.variables[1:nx]]

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
        if typeof(soln) <: ChebyshevSolutionStoch
            if soln.node_generator == chebyshev_nodes
                w = chebyshev_weights_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extrema
                w = chebyshev_weights_extrema_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extended
                w = chebyshev_weights_extended_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            end

            @sync @qthreads for t = 1:threads
                for j = t:threads:N
                    sub = ind2sub(j,Tuple(length.(grid)))

                    state = Array{T,1}(undef,nx)
                    for l = 1:nx
                        state[l] = grid[l][sub[l]]
                    end

                    variables[i][j] = chebyshev_evaluate(w,state,soln.order,domain)

                end
            end
        elseif typeof(soln) <: SmolyakSolutionStoch
            w = smolyak_weights_threaded(soln_variables[i],soln.grid,soln.multi_index,soln.domain)

            @sync @qthreads for t = 1:threads
                for j = t:threads:N
                    sub = ind2sub(j,Tuple(length.(grid)))

                    state = Array{T,1}(undef,nx)
                    for l = 1:nx
                        state[l] = grid[l][sub[l]]
                    end

                    variables[i][j] = smolyak_evaluate(w,state,soln.multi_index,domain)

                end
            end

        elseif typeof(soln) <: PiecewiseLinearSolutionStoch

            @sync @qthreads for t = 1:threads
                for j = t:threads:N
                    sub = ind2sub(j,Tuple(length.(grid)))

                    state = Array{T,1}(undef,nx)
                    for l = 1:nx
                        state[l] = grid[l][sub[l]]
                    end

                    variables[i][j] = piecewise_linear_evaluate(soln_variables[i],soln.nodes,state)

                end
            end
        end
    end

    weights        = Array{Array{T,nx},1}(undef,length(jumps_approximated))
    scaled_weights = Array{Array{T,nx},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        if node_generator == chebyshev_nodes
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extrema
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extrema(variables[jumps_approximated[i]],grid,order,domain)
            end
        elseif node_generator == chebyshev_extended
            for i = 1:length(jumps_approximated)
                weights[i] = chebyshev_weights_extended(variables[jumps_approximated[i]],grid,order,domain)
            end
        end

        for i = 1:length(jumps_approximated)
            for j = 1:ns
                index = [1:ndims(weights[i]);]
                index[1],index[j] = index[j],index[1]
                scaled_weights[i] = permutedims(integrals[j].*permutedims(weights[i],index),index)
            end
        end

        @sync @qthreads for t = 1:threads
            for i = t:threads:N
                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init  = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function(state,scaled_weights,order,domain,chebyshev_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end

            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionStoch([variables[ny+1:end];variables[1:ny]],weights,grid,order,domain,Matrix(k*k'),iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::SmolyakSchemeDet)

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    layer          = scheme.layer
    domain         = scheme.domain

    T = typeof(scheme.tol_fix_point_solver)

    grid, multi_ind = smolyak_grid(node_generator,nx,layer,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],N)
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],N)
    end

    weights = Array{Array{T,1},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(N)
    end

    init  = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:length(jumps_approximated)
            weights[i] = smolyak_weights(variables[jumps_approximated[i]],grid,multi_ind,domain)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function(grid[i,:],weights,multi_ind,domain,smolyak_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionDet([variables[ny+1:end];variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::SmolyakSchemeDet,threads::S) where {S <: Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    layer          = scheme.layer
    domain         = scheme.domain

    T = typeof(scheme.tol_fix_point_solver)

    grid, multi_ind = smolyak_grid(node_generator,nx,layer,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],N)
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],N)
    end

    weights = Array{Array{T,1},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(N)
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:length(jumps_approximated)
            weights[i] = smolyak_weights_threaded(variables[jumps_approximated[i]],grid,multi_ind,domain)
        end

        @sync @qthreads for t = 1:threads
            for i = t:threads:N

                init  = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function(grid[i,:],weights,multi_ind,domain,smolyak_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionDet([variables[ny+1:end];variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::SmolyakSchemeStoch)

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer          = scheme.layer
    domain         = scheme.domain

    T = typeof(scheme.tol_fix_point_solver)

    d = compute_linearization(model,initial_guess)
    k = -d[1:ns,2*nv+1:end]
    RHO = -d[1:ns,1:ns]
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    grid, multi_ind = smolyak_grid(node_generator,nx,layer,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],N)
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],N)
    end

    weights        = Array{Array{T,1},1}(undef,length(jumps_approximated))
    scaled_weights = Array{Array{T,1},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(N)
    end

    init  = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:length(jumps_approximated)
            weights[i]        = smolyak_weights(variables[jumps_approximated[i]],grid,multi_ind,domain)
            scaled_weights[i] = scale_weights(weights[i],weight_scale_factor)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function(grid[i,:],scaled_weights,multi_ind,domain,smolyak_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionStoch([variables[ny+1:end];variables[1:ny]],weights,grid,multi_ind,layer,domain,Matrix(k*k'),iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::SmolyakSchemeStoch,threads::S) where {S <: Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer          = scheme.layer
    domain         = scheme.domain

    T = typeof(scheme.tol_fix_point_solver)

    d = compute_linearization(model,initial_guess)
    k = -d[1:ns,2*nv+1:end]
    RHO = -d[1:ns,1:ns]
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    grid, multi_ind = smolyak_grid(node_generator,nx,layer,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],N)
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],N)
    end

    weights        = Array{Array{T,1},1}(undef,length(jumps_approximated))
    scaled_weights = Array{Array{T,1},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(N)
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:length(jumps_approximated)
            weights[i]        = smolyak_weights_threaded(variables[jumps_approximated[i]],grid,multi_ind,domain)
            scaled_weights[i] = scale_weights(weights[i],weight_scale_factor)
        end

        @sync @qthreads for t = 1:threads
            for i = t:threads:N

                init  = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function(grid[i,:],scaled_weights,multi_ind,domain,smolyak_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionStoch([variables[ny+1:end];variables[1:ny]],weights,grid,multi_ind,layer,domain,Matrix(k*k'),iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::SmolyakSchemeDet) where {R <: PerturbationSolutionDet}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    layer          = scheme.layer
    domain         = scheme.domain

    hbar = soln.hbar
    hx   = soln.hx
    gbar = soln.gbar
    gx   = soln.gx

    T = typeof(scheme.tol_fix_point_solver)

    grid, multi_ind = smolyak_grid(node_generator,nx,layer,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    for j = 1:N

        for i = 1:ny
            variables[i][j] = gbar[i] + (gx[i:i,:]*(grid[j,:] - hbar))[1]
        end
        for i = 1:nx
            variables[ny+i][j] = hbar[i] + (hx[i:i,:]*(grid[j,:] - hbar))[1]
        end

        if typeof(soln) <: SecondOrderSolutionDet
            for i = 1:ny
                variables[i][j] += (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))
            end
        end

        if typeof(soln) <: ThirdOrderSolutionDet
            for i = 1:ny
                variables[i][j] += (1/2)*((soln.gxx[i:i,:])*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((grid[j,:] - hbar),(grid[j,:] - hbar)),(grid[j,:] - hbar)))[1]
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*((soln.hxx[i:i,:])*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((grid[j,:] - hbar),(grid[j,:] - hbar)),(grid[j,:] - hbar)))[1]
            end
        end

    end

    weights = Array{Array{T,1},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(N)
    end

    init  = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:length(jumps_approximated)
            weights[i] = smolyak_weights(variables[jumps_approximated[i]],grid,multi_ind,domain)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function(grid[i,:],weights,multi_ind,domain,smolyak_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionDet([variables[ny+1:end];variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::SmolyakSchemeDet,threads) where {R <: PerturbationSolutionDet, S <: Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    layer          = scheme.layer
    domain         = scheme.domain

    hbar = soln.hbar
    hx   = soln.hx
    gbar = soln.gbar
    gx   = soln.gx

    T = typeof(scheme.tol_fix_point_solver)

    grid, multi_ind = smolyak_grid(node_generator,nx,layer,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    @sync @qthreads for t = 1:threads
        for j = t:threads:N

            for i = 1:ny
                variables[i][j] = gbar[i] + (gx[i:i,:]*(grid[j,:] - hbar))[1]
            end
            for i = 1:nx
                variables[ny+i][j] = hbar[i] + (hx[i:i,:]*(grid[j,:] - hbar))[1]
            end

            if typeof(soln) <: SecondOrderSolutionDet
                for i = 1:ny
                    variables[i][j] += (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))
                end
                for i = 1:nx
                    variables[ny+i][j] += (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))
                end
            end

            if typeof(soln) <: ThirdOrderSolutionDet
                for i = 1:ny
                    variables[i][j] += (1/2)*((soln.gxx[i:i,:])*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((grid[j,:] - hbar),(grid[j,:] - hbar)),(grid[j,:] - hbar)))[1]
                end
                for i = 1:nx
                    variables[ny+i][j] += (1/2)*((soln.hxx[i:i,:])*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((grid[j,:] - hbar),(grid[j,:] - hbar)),(grid[j,:] - hbar)))[1]
                end
            end
        end

    end

    weights = Array{Array{T,1},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(N)
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:length(jumps_approximated)
            weights[i] = smolyak_weights_threaded(variables[jumps_approximated[i]],grid,multi_ind,domain)
        end

        @sync @qthreads for t = 1:threads
            for i = t:threads:N

                init  = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function(grid[i,:],weights,multi_ind,domain,smolyak_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionDet([variables[ny+1:end];variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::SmolyakSchemeStoch) where {R <: PerturbationSolutionStoch}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer          = scheme.layer
    domain         = scheme.domain

    hbar = soln.hbar
    hx   = soln.hx
    k    = soln.k
    gbar = soln.gbar
    gx   = soln.gx

    RHO = hx[1:ns,1:ns]

    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        if typeof(soln) <: FirstOrderSolutionStoch
            soln_fo = soln
        else
            soln_fo = FirstOrderSolutionStoch(hbar,hx,k,gbar,gx,soln.sigma,soln.grc,soln.soln_type)
        end

        state_vars, jump_vars = compute_variances(soln_fo)

        domain = Matrix([hbar + 3*sqrt.(diag(state_vars)) hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
    end

    grid, multi_ind = smolyak_grid(node_generator,nx,layer,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    for j = 1:N

        for i = 1:ny
            variables[i][j] = gbar[i] + (gx[i:i,:]*(grid[j,:] - hbar))[1]
        end
        for i = 1:nx
            variables[ny+i][j] = hbar[i] + (hx[i:i,:]*(grid[j,:] - hbar))[1]
        end

        if typeof(soln) <: SecondOrderSolutionStoch
            for i = 1:ny
                variables[i][j] += (1/2)*soln.gss[i] + (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*soln.hss[i] + (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))
            end
        end

        if typeof(soln) <: ThirdOrderSolutionStoch
            for i = 1:ny
                variables[i][j] += (1/2)*soln.gss[i] + (1/2)*((soln.gxx[i:i,:])*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))[1] + (3/6)*(soln.gssx[i:i,:]*(grid[j,:] - hbar))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((grid[j,:] - hbar),(grid[j,:] - hbar)),(grid[j,:] - hbar)))[1]
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*soln.hss[i] + (1/2)*((soln.hxx[i:i,:])*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))[1] + (3/6)*(soln.hssx[i:i,:]*(grid[j,:] - hbar))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((grid[j,:] - hbar),(grid[j,:] - hbar)),(grid[j,:] - hbar)))[1]
            end
        end

    end

    weights        = Array{Array{T,1},1}(undef,length(jumps_approximated))
    scaled_weights = Array{Array{T,1},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(N)
    end

    init  = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:length(jumps_approximated)
            weights[i]        = smolyak_weights(variables[jumps_approximated[i]],grid,multi_ind,domain)
            scaled_weights[i] = scale_weights(weights[i],weight_scale_factor)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function(grid[i,:],scaled_weights,multi_ind,domain,smolyak_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionStoch([variables[ny+1:end];variables[1:ny]],weights,grid,multi_ind,layer,domain,Matrix(k[1:ns,:]*k[1:ns,:]'),iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::SmolyakSchemeStoch,threads::S) where {R <: PerturbationSolutionStoch, S <: Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer          = scheme.layer
    domain         = scheme.domain

    hbar = soln.hbar
    hx   = soln.hx
    k    = soln.k
    gbar = soln.gbar
    gx   = soln.gx

    RHO = hx[1:ns,1:ns]

    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        if typeof(soln) <: FirstOrderSolutionStoch
            soln_fo = soln
        else
            soln_fo = FirstOrderSolutionStoch(hbar,hx,k,gbar,gx,soln.sigma,soln.grc,soln.soln_type)
        end

        state_vars, jump_vars = compute_variances(soln_fo)

        domain = Matrix([hbar + 3*sqrt.(diag(state_vars)) hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
    end

    grid, multi_ind = smolyak_grid(node_generator,nx,layer,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    @sync @qthreads for t = 1:threads
        for j = t:threads:N

            for i = 1:ny
                variables[i][j] = gbar[i] + (gx[i:i,:]*(grid[j,:] - hbar))[1]
            end
            for i = 1:nx
                variables[ny+i][j] = hbar[i] + (hx[i:i,:]*(grid[j,:] - hbar))[1]
            end

            if typeof(soln) <: SecondOrderSolutionStoch
                for i = 1:ny
                    variables[i][j] += (1/2)*soln.gss[i] + (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))
                end
                for i = 1:nx
                    variables[ny+i][j] += (1/2)*soln.hss[i] + (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))
                end
            end

            if typeof(soln) <: ThirdOrderSolutionStoch
                for i = 1:ny
                    variables[i][j] += (1/2)*soln.gss[i] + (1/2)*((soln.gxx[i:i,:])*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))[1] + (3/6)*(soln.gssx[i:i,:]*(grid[j,:] - hbar))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((grid[j,:] - hbar),(grid[j,:] - hbar)),(grid[j,:] - hbar)))[1]
                end
                for i = 1:nx
                    variables[ny+i][j] += (1/2)*soln.hss[i] + (1/2)*((soln.hxx[i:i,:])*kron((grid[j,:] - hbar),(grid[j,:] - hbar)))[1] + (3/6)*(soln.hssx[i:i,:]*(grid[j,:] - hbar))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((grid[j,:] - hbar),(grid[j,:] - hbar)),(grid[j,:] - hbar)))[1]
                end
            end
        end
    end

    weights        = Array{Array{T,1},1}(undef,length(jumps_approximated))
    scaled_weights = Array{Array{T,1},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(N)
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:length(jumps_approximated)
            weights[i]        = smolyak_weights_threaded(variables[jumps_approximated[i]],grid,multi_ind,domain)
            scaled_weights[i] = scale_weights(weights[i],weight_scale_factor)
        end

        @sync @qthreads for t = 1:threads
            for i = t:threads:N

                init  = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function(grid[i,:],scaled_weights,multi_ind,domain,smolyak_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionStoch([variables[ny+1:end];variables[1:ny]],weights,grid,multi_ind,layer,domain,Matrix(k[1:ns,:]*k[1:ns,:]'),iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::SmolyakSchemeDet) where {R <: ProjectionSolutionDet}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    layer          = scheme.layer
    domain         = scheme.domain

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        domain = soln.domain
    end

    grid, multi_ind = smolyak_grid(node_generator,nx,layer,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    soln_variables = [soln.variables[nx+1:end];soln.variables[1:nx]]

    for i = 1:nv
        variables[i] = zeros(N)
        if typeof(soln) <: ChebyshevSolutionDet

            if soln.node_generator == chebyshev_nodes
                w = chebyshev_weights_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extrema
                w = chebyshev_weights_extrema_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extended
                w = chebyshev_weights_extended_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            end

            for j = 1:N

                variables[i][j] = chebyshev_evaluate(w,grid[j,:],soln.order,domain)

            end

        elseif typeof(soln) <: SmolyakSolutionDet

            w = smolyak_weights(soln_variables[i],soln.grid,soln.multi_index,soln.domain)

            for j = 1:N

                variables[i][j] = smolyak_evaluate(w,grid[j,:],soln.multi_index,soln.domain)

            end

        elseif typeof(soln) <: PiecewiseLinearSolutionDet

            for j = 1:N

                variables[i][j] = piecewise_linear_evaluate(soln_variables[i],soln.nodes,grid[j,:])

            end

        end
    end

    weights = Array{Array{T,1},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(N)
    end

    init  = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:length(jumps_approximated)
            weights[i] = smolyak_weights(variables[jumps_approximated[i]],grid,multi_ind,domain)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function(grid[i,:],weights,multi_ind,domain,smolyak_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionDet([variables[ny+1:end];variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::SmolyakSchemeDet,threads::S) where {R <: ProjectionSolutionDet, S <: Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    layer          = scheme.layer
    domain         = scheme.domain

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        domain = soln.domain
    end

    grid, multi_ind = smolyak_grid(node_generator,nx,layer,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    soln_variables = [soln.variables[nx+1:end];soln.variables[1:nx]]

    for i = 1:nv
        variables[i] = zeros(N)
        if typeof(soln) <: ChebyshevSolutionDet

            if soln.node_generator == chebyshev_nodes
                w = chebyshev_weights_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extrema
                w = chebyshev_weights_extrema_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extended
                w = chebyshev_weights_extended_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            end

            @sync @qthreads for t = 1:threads
                for j = t:threads:N

                    variables[i][j] = chebyshev_evaluate(w,grid[j,:],soln.order,domain)

                end
            end

        elseif typeof(soln) <: SmolyakSolutionDet

            w = smolyak_weights_threaded(soln_variables[i],soln.grid,soln.multi_index,soln.domain)

            @sync @qthreads for t = 1:threads
                for j = t:threads:N

                    variables[i][j] = smolyak_evaluate(w,grid[j,:],soln.multi_index,soln.domain)

                end
            end

        elseif typeof(soln) <: PiecewiseLinearSolutionDet

            @sync @qthreads for t = 1:threads
                for j = t:threads:N

                    variables[i][j] = piecewise_linear_evaluate(soln_variables[i],soln.nodes,grid[j,:])

                end
            end

        end
    end

    weights = Array{Array{T,1},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(N)
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:length(jumps_approximated)
            weights[i] = smolyak_weights_threaded(variables[jumps_approximated[i]],grid,multi_ind,domain)
        end

        @sync @qthreads for t = 1:threads
            for i = t:threads:N

                init  = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function(grid[i,:],weights,multi_ind,domain,smolyak_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionDet([variables[ny+1:end];variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::SmolyakSchemeStoch) where {R <: ProjectionSolutionStoch}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer          = scheme.layer
    domain         = scheme.domain

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        domain = soln.domain
    end

    d = compute_linearization(model,initial_guess)
    k = -d[1:ns,2*nv+1:end]
    RHO = -d[1:ns,1:ns]
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    grid, multi_ind = smolyak_grid(node_generator,nx,layer,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    soln_variables = [soln.variables[nx+1:end];soln.variables[1:nx]]

    for i = 1:nv
        variables[i] = zeros(N)
        if typeof(soln) <: ChebyshevSolutionStoch

            if soln.node_generator == chebyshev_nodes
                w = chebyshev_weights_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extrema
                w = chebyshev_weights_extrema_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extended
                w = chebyshev_weights_extended_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            end

            for j = 1:N

                variables[i][j] = chebyshev_evaluate(w,grid[j,:],soln.order,domain)

            end

        elseif typeof(soln) <: SmolyakSolutionStoch

            w = smolyak_weights(soln_variables[i],soln.grid,soln.multi_index,soln.domain)

            for j = 1:N

                variables[i][j] = smolyak_evaluate(w,grid[j,:],soln.multi_index,soln.domain)

            end

        elseif typeof(soln) <: PiecewiseLinearSolutionStoch

            for j = 1:N

                variables[i][j] = piecewise_linear_evaluate(soln_variables[i],soln.nodes,grid[j,:])

            end

        end
    end

    weights        = Array{Array{T,1},1}(undef,length(jumps_approximated))
    scaled_weights = Array{Array{T,1},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(N)
    end

    init  = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:length(jumps_approximated)
            weights[i]        = smolyak_weights(variables[jumps_approximated[i]],grid,multi_ind,domain)
            scaled_weights[i] = scale_weights(weights[i],weight_scale_factor)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function(grid[i,:],scaled_weights,multi_ind,domain,smolyak_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionStoch([variables[ny+1:end];variables[1:ny]],weights,grid,multi_ind,layer,domain,Matrix(k*k'),iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::SmolyakSchemeStoch,threads::S) where {R <: ProjectionSolutionStoch, S <: Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer          = scheme.layer
    domain         = scheme.domain

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        domain = soln.domain
    end

    d = compute_linearization(model,initial_guess)
    k = -d[1:ns,2*nv+1:end]
    RHO = -d[1:ns,1:ns]
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    grid, multi_ind = smolyak_grid(node_generator,nx,layer,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    soln_variables = [soln.variables[nx+1:end];soln.variables[1:nx]]

    for i = 1:nv
        variables[i] = zeros(N)
        if typeof(soln) <: ChebyshevSolutionStoch

            if soln.node_generator == chebyshev_nodes
                w = chebyshev_weights_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extrema
                w = chebyshev_weights_extrema_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extended
                w = chebyshev_weights_extended_threaded(soln_variables[i],soln.nodes,soln.order,soln.domain)
            end

            @sync @qthreads for t = 1:threads
                for j = t:threads:N

                    variables[i][j] = chebyshev_evaluate(w,grid[j,:],soln.order,domain)

                end
            end

        elseif typeof(soln) <: SmolyakSolutionStoch

            w = smolyak_weights_threaded(soln_variables[i],soln.grid,soln.multi_index,soln.domain)

            @sync @qthreads for t = 1:threads
                for j = t:threads:N

                    variables[i][j] = smolyak_evaluate(w,grid[j,:],soln.multi_index,soln.domain)

                end
            end

        elseif typeof(soln) <: PiecewiseLinearSolutionStoch

            @sync @qthreads for t = 1:threads
                for j = t:threads:N

                    variables[i][j] = piecewise_linear_evaluate(soln_variables[i],soln.nodes,grid[j,:])

                end
            end

        end
    end

    weights        = Array{Array{T,1},1}(undef,length(jumps_approximated))
    scaled_weights = Array{Array{T,1},1}(undef,length(jumps_approximated))

    new_variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(N)
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:length(jumps_approximated)
            weights[i]        = smolyak_weights_threaded(variables[jumps_approximated[i]],grid,multi_ind,domain)
            scaled_weights[i] = scale_weights(weights[i],weight_scale_factor)
        end

        @sync @qthreads for t = 1:threads
            for i = t:threads:N

                init  = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function(grid[i,:],scaled_weights,multi_ind,domain,smolyak_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionStoch([variables[ny+1:end];variables[1:ny]],weights,grid,multi_ind,layer,domain,Matrix(k*k'),iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::PiecewiseLinearSchemeDet)

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_number    = scheme.node_number
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    state = Array{T,1}(undef,nx)

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],Tuple(length.(grid)))
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],Tuple(length.(grid)))
    end

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    init = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:N
            sub = ind2sub(i,Tuple(length.(grid)))

            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_piecewise(variables,grid,state,piecewise_linear_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionDet([variables[ny+1:end];variables[1:ny]],grid,domain,iters)

    return soln

end

function solve_nonlinear(model::REModel,scheme::PiecewiseLinearSchemeDet,threads::S) where {S <: Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_number    = scheme.node_number
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],Tuple(length.(grid)))
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],Tuple(length.(grid)))
    end

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        @sync @qthreads for t = 1:threads
            for i = t:threads:N
                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init  = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_piecewise(variables,grid,state,piecewise_linear_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionDet([variables[ny+1:end];variables[1:ny]],grid,domain,iters)

    return soln

end

function solve_nonlinear(model::REModel,scheme::PiecewiseLinearSchemeStoch)

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_number    = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    d = compute_linearization(model,initial_guess)
    k = -d[1:ns,2*nv+1:end]
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = ones(nx)
    for i = 1:ns
        integrals[i] = compute_piecewise_linear_integrals(eps_nodes,eps_weights,k[i,i])
    end

    state = Array{T,1}(undef,nx)

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],Tuple(length.(grid)))
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],Tuple(length.(grid)))
    end

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    init  = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:N
            sub = ind2sub(i,Tuple(length.(grid)))

            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_piecewise(variables,grid,state,integrals,piecewise_linear_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionStoch([variables[ny+1:end];variables[1:ny]],grid,domain,Matrix(k*k'),iters)

    return soln

end

function solve_nonlinear(model::REModel,scheme::PiecewiseLinearSchemeStoch,threads::S) where {S <: Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess  = scheme.initial_guess
    node_number    = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    d = compute_linearization(model,initial_guess)
    k = -d[1:ns,2*nv+1:end]
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = ones(nx)
    for i = 1:ns
        integrals[i] = compute_piecewise_linear_integrals(eps_nodes,eps_weights,k[i,i])
    end

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],Tuple(length.(grid)))
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],Tuple(length.(grid)))
    end

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        @sync @qthreads for t = 1:threads
            for i = t:threads:N
                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init  = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_piecewise(variables,grid,state,integrals,piecewise_linear_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionStoch([variables[ny+1:end];variables[1:ny]],grid,domain,Matrix(k*k'),iters)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::PiecewiseLinearSchemeDet) where {R <: PerturbationSolutionDet}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    initial_guess  = scheme.initial_guess
    node_number    = scheme.node_number
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    hbar = soln.hbar
    hx   = soln.hx
    gbar = soln.gbar
    gx   = soln.gx

    T = typeof(scheme.tol_fix_point_solver)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    state = Array{T,1}(undef,nx)

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    for j = 1:N
        sub = ind2sub(j,Tuple(length.(grid)))

        for i = 1:nx
            state[i] = grid[i][sub[i]]
        end

        for i = 1:ny
            variables[i][j] = gbar[i] + (gx[i:i,:]*(state - hbar))[1]
        end
        for i = 1:nx
            variables[ny+i][j] = hbar[i] + (hx[i:i,:]*(state - hbar))[1]
        end

        if typeof(soln) <: SecondOrderSolutionDet
            for i = 1:ny
                variables[i][j] += (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
            end
        end

        if typeof(soln) <: ThirdOrderSolutionDet
            for i = 1:ny
                variables[i][j] += (1/2)*((soln.gxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*((soln.hxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
            end
        end

    end

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    init  = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:N
            sub = ind2sub(i,Tuple(length.(grid)))

            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_piecewise(variables,grid,state,piecewise_linear_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionDet([variables[ny+1:end];variables[1:ny]],grid,domain,iters)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::PiecewiseLinearSchemeDet,threads::S) where {R <: PerturbationSolutionDet, S <: Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    initial_guess  = scheme.initial_guess
    node_number    = scheme.node_number
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    hbar = soln.hbar
    hx   = soln.hx
    gbar = soln.gbar
    gx   = soln.gx

    T = typeof(scheme.tol_fix_point_solver)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    @sync @qthreads for t = 1:threads
        for j = t:threads:N
            sub = ind2sub(j,Tuple(length.(grid)))

            state = Array{T,1}(undef,nx)
            for i = 1:nx
                state[i] = grid[i][sub[i]]
            end

            for i = 1:ny
                variables[i][j] = gbar[i] + (gx[i:i,:]*(state - hbar))[1]
            end
            for i = 1:nx
                variables[ny+i][j] = hbar[i] + (hx[i:i,:]*(state - hbar))[1]
            end

            if typeof(soln) <: SecondOrderSolutionDet
                for i = 1:ny
                    variables[i][j] += (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
                end
                for i = 1:nx
                    variables[ny+i][j] += (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
                end
            end

            if typeof(soln) <: ThirdOrderSolutionDet
                for i = 1:ny
                    variables[i][j] += (1/2)*((soln.gxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
                end
                for i = 1:nx
                    variables[ny+i][j] += (1/2)*((soln.hxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
                end
            end
        end
    end

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        @sync @qthreads for t = 1:threads
            for i = t:threads:N
                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init  = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_piecewise(variables,grid,state,piecewise_linear_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionDet([variables[ny+1:end];variables[1:ny]],grid,domain,iters)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::PiecewiseLinearSchemeStoch) where {R <: PerturbationSolutionStoch}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    initial_guess  = scheme.initial_guess
    node_number    = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    hbar = soln.hbar
    hx   = soln.hx
    k    = soln.k
    gbar = soln.gbar
    gx   = soln.gx

    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        if typeof(soln) <: FirstOrderSolutionStoch
            soln_fo = soln
        else
            soln_fo = FirstOrderSolutionStoch(hbar,hx,k,gbar,gx,soln.sigma,soln.grc,soln.soln_type)
        end

        state_vars, jump_vars = compute_variances(soln_fo)

        domain = Matrix([hbar + 3*sqrt.(diag(state_vars)) hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = ones(nx)
    for i = 1:ns
        integrals[i] = compute_piecewise_linear_integrals(eps_nodes,eps_weights,k[i,i])
    end

    state = Array{T,1}(undef,nx)

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    for j = 1:N
        sub = ind2sub(j,Tuple(length.(grid)))

        for i = 1:nx
            state[i] = grid[i][sub[i]]
        end

        for i = 1:ny
            variables[i][j] = gbar[i] + (gx[i:i,:]*(state - hbar))[1]
        end
        for i = 1:nx
            variables[ny+i][j] = hbar[i] + (hx[i:i,:]*(state - hbar))[1]
        end

        if typeof(soln) <: SecondOrderSolutionStoch
            for i = 1:ny
                variables[i][j] += (1/2)*soln.gss[i] + (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*soln.hss[i] + (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
            end
        end

        if typeof(soln) <: ThirdOrderSolutionStoch
            for i = 1:ny
                variables[i][j] += (1/2)*soln.gss[i] + (1/2)*((soln.gxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (3/6)*(soln.gssx[i:i,:]*(state - hbar))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
            end
            for i = 1:nx
                variables[ny+i][j] += (1/2)*soln.hss[i] + (1/2)*((soln.hxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (3/6)*(soln.hssx[i:i,:]*(state - hbar))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
            end
        end

    end

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    init  = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:N
            sub = ind2sub(i,Tuple(length.(grid)))

            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_piecewise(variables,grid,state,integrals,piecewise_linear_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionStoch([variables[ny+1:end];variables[1:ny]],grid,domain,Matrix(k[1:ns,:]*k[1:ns,:]'),iters)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::PiecewiseLinearSchemeStoch,threads::S) where {R <: PerturbationSolutionStoch, S <: Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    initial_guess  = scheme.initial_guess
    node_number    = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    hbar = soln.hbar
    hx   = soln.hx
    k    = soln.k
    gbar = soln.gbar
    gx   = soln.gx

    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        if typeof(soln) <: FirstOrderSolutionStoch
            soln_fo = soln
        else
            soln_fo = FirstOrderSolutionStoch(hbar,hx,k,gbar,gx,soln.sigma,soln.grc,soln.soln_type)
        end

        state_vars, jump_vars = compute_variances(soln_fo)

        domain = Matrix([hbar + 3*sqrt.(diag(state_vars)) hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = ones(nx)
    for i = 1:ns
        integrals[i] = compute_piecewise_linear_integrals(eps_nodes,eps_weights,k[i,i])
    end

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    @sync @qthreads for t = 1:threads
        for j = t:threads:N
            sub = ind2sub(j,Tuple(length.(grid)))

            state = Array{T,1}(undef,nx)
            for i = 1:nx
                state[i] = grid[i][sub[i]]
            end

            for i = 1:ny
                variables[i][j] = gbar[i] + (gx[i:i,:]*(state - hbar))[1]
            end
            for i = 1:nx
                variables[ny+i][j] = hbar[i] + (hx[i:i,:]*(state - hbar))[1]
            end

            if typeof(soln) <: SecondOrderSolutionStoch
                for i = 1:ny
                    variables[i][j] += (1/2)*soln.gss[i] + (1/2)*sum(vec(soln.gxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
                end
                for i = 1:nx
                    variables[ny+i][j] += (1/2)*soln.hss[i] + (1/2)*sum(vec(soln.hxx[(i-1)*nx+1:i*nx,:]).*kron((state - hbar),(state - hbar)))
                end
            end

            if typeof(soln) <: ThirdOrderSolutionStoch
                for i = 1:ny
                    variables[i][j] += (1/2)*soln.gss[i] + (1/2)*((soln.gxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (3/6)*(soln.gssx[i:i,:]*(state - hbar))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
                end
                for i = 1:nx
                    variables[ny+i][j] += (1/2)*soln.hss[i] + (1/2)*((soln.hxx[i:i,:])*kron((state - hbar),(state - hbar)))[1] + (3/6)*(soln.hssx[i:i,:]*(state - hbar))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state - hbar),(state - hbar)),(state-hbar)))[1]
                end
            end
        end
    end

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        @sync @qthreads for t = 1:threads
            for i = t:threads:N
                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init  = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_piecewise(variables,grid,state,integrals,piecewise_linear_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionStoch([variables[ny+1:end];variables[1:ny]],grid,domain,Matrix(k[1:ns,:]*k[1:ns,:]'),iters)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::PiecewiseLinearSchemeDet) where {R <: ProjectionSolutionDet}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    initial_guess  = scheme.initial_guess
    node_number    = scheme.node_number
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        domain = soln.domain
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    state = Array{T,1}(undef,nx)

    N = prod(length.(grid))

    soln_variables = [soln.variables[nx+1:end];soln.variables[1:nx]]

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
        if typeof(soln) <: ChebyshevSolutionDet

            if soln.node_generator == chebyshev_nodes
                w = chebyshev_weights(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extrema
                w = chebyshev_weights_extrema(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extended
                w = chebyshev_weights_extended(soln_variables[i],soln.nodes,soln.order,soln.domain)
            end

            for j = 1:N
                sub = ind2sub(j,Tuple(length.(grid)))

                for l = 1:nx
                    state[l] = grid[l][sub[l]]
                end

                variables[i][j] = chebyshev_evaluate(w,state,soln.order,domain)

            end
        elseif typeof(soln) <: SmolyakSolutionDet
            w = smolyak_weights(soln_variables[i],soln.grid,soln.multi_index,soln.domain)

            for j = 1:N
                sub = ind2sub(j,Tuple(length.(grid)))

                for l = 1:nx
                    state[l] = grid[l][sub[l]]
                end

                variables[i][j] = smolyak_evaluate(w,state,soln.multi_index,domain)

            end

        elseif typeof(soln) <: PiecewiseLinearSolutionDet

            variables[i] = PiecewiseLinearApprox.grid_reshape(soln_variables[i],tuple(grid...))

        end
    end

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    init  = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:N
            sub = ind2sub(i,Tuple(length.(grid)))

            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_piecewise(variables,grid,state,piecewise_linear_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionDet([variables[ny+1:end];variables[1:ny]],grid,domain,iters)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::PiecewiseLinearSchemeDet,threads::S) where {R <: ProjectionSolutionDet, S <: Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    initial_guess  = scheme.initial_guess
    node_number    = scheme.node_number
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        domain = soln.domain
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    N = prod(length.(grid))

    soln_variables = [soln.variables[nx+1:end];soln.variables[1:nx]]

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
        if typeof(soln) <: ChebyshevSolutionDet

            if soln.node_generator == chebyshev_nodes
                w = chebyshev_weights(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extrema
                w = chebyshev_weights_extrema(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extended
                w = chebyshev_weights_extended(soln_variables[i],soln.nodes,soln.order,soln.domain)
            end

            @sync @qthreads for t = 1:threads
                for j = t:threads:N
                    sub = ind2sub(j,Tuple(length.(grid)))

                    state = Array{T,1}(undef,nx)
                    for l = 1:nx
                        state[l] = grid[l][sub[l]]
                    end

                    variables[i][j] = chebyshev_evaluate(w,state,soln.order,domain)

                end
            end

        elseif typeof(soln) <: SmolyakSolutionDet

            w = smolyak_weights_threaded(soln_variables[i],soln.grid,soln.multi_index,soln.domain)

            @sync @qthreads for t = 1:threads
                for j = t:threads:N
                    sub = ind2sub(j,Tuple(length.(grid)))

                    state = Array{T,1}(undef,nx)
                    for l = 1:nx
                        state[l] = grid[l][sub[l]]
                    end

                    variables[i][j] = smolyak_evaluate(w,state,soln.multi_index,domain)

                end
            end

        elseif typeof(soln) <: PiecewiseLinearSolutionDet

            variables[i] = PiecewiseLinearApprox.grid_reshape(soln_variables[i],tuple(grid...))

        end
    end

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        @sync @qthreads for t = 1:threads
            for i = t:threads:N
                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init  = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_piecewise(variables,grid,state,piecewise_linear_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionDet([variables[ny+1:end];variables[1:ny]],grid,domain,iters)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::PiecewiseLinearSchemeStoch) where {R <: ProjectionSolutionStoch}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    initial_guess  = scheme.initial_guess
    node_number    = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        domain = soln.domain
    end

    d = compute_linearization(model,initial_guess)
    k = -d[1:ns,2*nv+1:end]
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = ones(nx)
    for i = 1:ns
        integrals[i] = compute_piecewise_linear_integrals(eps_nodes,eps_weights,k[i,i])
    end

    state = Array{T,1}(undef,nx)

    N = prod(length.(grid))

    soln_variables = [soln.variables[nx+1:end];soln.variables[1:nx]]

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
        if typeof(soln) <: ChebyshevSolutionStoch

            if soln.node_generator == chebyshev_nodes
                w = chebyshev_weights(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extrema
                w = chebyshev_weights_extrema(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extended
                w = chebyshev_weights_extended(soln_variables[i],soln.nodes,soln.order,soln.domain)
            end

            for j = 1:N
                sub = ind2sub(j,Tuple(length.(grid)))

                for l = 1:nx
                    state[l] = grid[l][sub[l]]
                end

                variables[i][j] = chebyshev_evaluate(w,state,soln.order,domain)

            end

        elseif typeof(soln) <: SmolyakSolutionStoch

            w = smolyak_weights(soln_variables[i],soln.grid,soln.multi_index,soln.domain)

            for j = 1:N
                sub = ind2sub(j,Tuple(length.(grid)))

                for l = 1:nx
                    state[l] = grid[l][sub[l]]
                end

                variables[i][j] = smolyak_evaluate(w,state,soln.multi_index,domain)

            end

        elseif typeof(soln) <: PiecewiseLinearSolutionStoch

            variables[i] = PiecewiseLinearApprox.grid_reshape(soln_variables[i],tuple(grid...))

        end
    end

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    init  = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        for i = 1:N
            sub = ind2sub(i,Tuple(length.(grid)))

            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_piecewise(variables,grid,state,integrals,piecewise_linear_evaluate)
            nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionStoch([variables[ny+1:end];variables[1:ny]],grid,domain,Matrix(k*k'),iters)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::PiecewiseLinearSchemeStoch,threads::S) where {R <: ProjectionSolutionStoch, S <: Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    initial_guess  = scheme.initial_guess
    node_number    = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    domain         = scheme.domain

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.tol_fix_point_solver)

    if domain == []
        domain = soln.domain
    end

    d = compute_linearization(model,initial_guess)
    k = -d[1:ns,2*nv+1:end]
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),scheme.tol_variables)) != ns-1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via projection methods")
        end
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = ones(nx)
    for i = 1:ns
        integrals[i] = compute_piecewise_linear_integrals(eps_nodes,eps_weights,k[i,i])
    end

    N = prod(length.(grid))

    soln_variables = [soln.variables[nx+1:end];soln.variables[1:nx]]

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
        if typeof(soln) <: ChebyshevSolutionStoch

            if soln.node_generator == chebyshev_nodes
                w = chebyshev_weights(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extrema
                w = chebyshev_weights_extrema(soln_variables[i],soln.nodes,soln.order,soln.domain)
            elseif soln.node_generator == chebyshev_extended
                w = chebyshev_weights_extended(soln_variables[i],soln.nodes,soln.order,soln.domain)
            end

            @sync @qthreads for t = 1:threads
                for j = t:threads:N
                    sub = ind2sub(j,Tuple(length.(grid)))

                    state = Array{T,1}(undef,nx)
                    for l = 1:nx
                        state[l] = grid[l][sub[l]]
                    end

                    variables[i][j] = chebyshev_evaluate(w,state,soln.order,domain)

                end
            end
        elseif typeof(soln) <: SmolyakSolutionStoch

            w = smolyak_weights_threaded(soln_variables[i],soln.grid,soln.multi_index,soln.domain)

            @sync @qthreads for t = 1:threads
                for j = t:threads:N
                    sub = ind2sub(j,Tuple(length.(grid)))

                    state = Array{T,1}(undef,nx)
                    for l = 1:nx
                        state[l] = grid[l][sub[l]]
                    end

                    variables[i][j] = smolyak_evaluate(w,state,soln.multi_index,domain)

                end
            end

        elseif typeof(soln) <: PiecewiseLinearSolutionStoch

            variables[i] = PiecewiseLinearApprox.grid_reshape(soln_variables[i],tuple(grid...))

        end
    end

    new_variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        new_variables[i] = zeros(Tuple(length.(grid)))
    end

    iters = 0
    len = Inf
    while len > scheme.tol_variables && iters <= scheme.maxiters

        @sync @qthreads for t = 1:threads
            for i = t:threads:N
                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init  = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_piecewise(variables,grid,state,integrals,piecewise_linear_evaluate)
                nlsoln = nlsolve(projection_equations, init, xtol = scheme.tol_fix_point_solver, iterations = scheme.maxiters, inplace = :true)
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = ThreadsX.maximum(abs,new_variables[1]-variables[1])
        for j = 2:nv
           len = max(ThreadsX.maximum(abs,new_variables[j]-variables[j]),len)
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionStoch([variables[ny+1:end];variables[1:ny]],grid,domain,Matrix(k*k'),iters)

    return soln

end

function solve_model(model::REModel,scheme::PerturbationScheme)

    if scheme.order == "first"
        soln = solve_first_order(model, scheme)
        return soln
    elseif scheme.order == "second"
        soln = solve_second_order(model, scheme)
        return soln
    elseif scheme.order == "third"
        soln = solve_third_order(model, scheme)
        return soln
    else
        error("The chosen order has not been implemented")
    end

end

function solve_model(model::REModel,scheme::P) where {P <: ProjectionScheme}

    soln = solve_nonlinear(model,scheme)
    return soln

end

function solve_model(model::REModel,scheme::P,threads::S) where {P <: ProjectionScheme, S <: Integer}

    if threads < 0
        error("Number of threads cannot be negative")
    elseif threads == 0
        soln = solve_nonlinear(model,scheme)
        return soln
    else
        soln = solve_nonlinear(model,scheme,threads)
        return soln
    end

end

function solve_model(model::REModel,soln::ModelSolution,scheme::P) where {P <: ProjectionScheme}

    soln = solve_nonlinear(model,soln,scheme)
    return soln

end

function solve_model(model::REModel,soln::ModelSolution,scheme::P,threads::S) where {P <: ProjectionScheme, S <: Integer}

    if threads < 0
        error("Number of threads cannot be negative")
    elseif threads == 0
        soln = solve_nonlinear(model,soln,scheme)
        return soln
    else
        soln = solve_nonlinear(model,soln,scheme,threads)
        return soln
    end

end
