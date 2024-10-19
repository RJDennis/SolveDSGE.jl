###################### Solution functions #############################

"""
Computes the model's deterministic steady state.

Exported function.
"""
function compute_steady_state(model::REModel,x::Array{T,1},tol::T,maxiters::S,method::Symbol = :newton) where {T<:Real,S<:Integer}
#function compute_steady_state(model::REModel,x::Array{T,1},xtol::T,ftol::T,maxiters::S,method::Symbol = :lm_ar) where {T<:Real,S<:Integer}

    if length(x) > model.number_equations
        error("The initialization has too many elements.")
    elseif length(x) < model.number_equations
        error("The initialization has too few elements.")
    end

    #equations = model.static_function
    #nlsoln = nlsolve(equations, x, xtol = tol, iterations = maxiters, autodiff = :forward, inplace = :false)
    equations = model.nlsolve_static_function
    nlsoln = nlsolve(equations,x,ftol = tol,xtol = tol,iterations = maxiters,autodiff = :forward,inplace = :true,method = method)
    #nlsoln = nlboxsolve(equations,x,xtol = xtol,ftol = ftol,iterations = maxiters,method = method)

    return nlsoln

end

"""
Computes the Jacobian of the model's equations using automatic derivatives.

Internal function; not exposed to users.
"""
function compute_linearization(model::REModel,steady_state::Array{T,1}) where {T<:Real}

    equations = model.dynamic_function
    ns = model.number_shocks

    x = [steady_state; steady_state; zeros(ns)]    
    d = ForwardDiff.jacobian(equations,x)

    return d

end

"""
Implements Klein's (2000) method to solve a model to first-order accuracy.

Exported function.
"""
function solve_first_order(model::REModel,scheme::PerturbationScheme)

    if scheme.order != "first"
        error("A first order perturbation must be specified.")
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

"""
Implements Klein's (2000) method to solve a deterministic model to first-order accuracy.

Internal function; not exposed to users.
"""
function solve_first_order_det(model::REModel,scheme::PerturbationScheme)

    nx = model.number_states
    ny = model.number_jumps
    nv = model.number_variables

    steady_state = scheme.steady_state
    cutoff = scheme.cutoff

    # Compute the first-order deterministic terms

    d = compute_linearization(model,steady_state)

    @views a = d[1:nv,1:nv]
    @views b = -d[1:nv,nv+1:2*nv]

    r = schur(complex(a),complex(b))

    # Reorder the generalized eigenvalues so that those with modulus greater
    # than "cutoff" reside at the bottom.

    sel = (abs.(diag(r.S)./diag(r.T)) .< cutoff)
    ordschur!(r,sel)

    s = r.S
    t = r.T
    #q = r.Q'  # So now q*a*z = s and q*b*z = t
    z = r.Z

    # Construct the rational expectations equilibrium by eliminating the
    # unstable dynamics

    @views z11 = z[1:nx,1:nx]
    @views z21 = z[nx+1:nv,1:nx]
    @views t11 = t[1:nx,1:nx]
    @views s11 = s[1:nx,1:nx]

    if rank(z11) != nx
        error("The model does not have a stable solution.")
    end

    hx = real((z11/t11)*(s11/z11))
    gx = real(z21/z11)

    # Calculate the number of unstable eigenvalues

    grc = nv - sum(sel)

    soln_type = "determinate"
    if grc < ny
        soln_type = "indeterminate"
    elseif grc > ny
        soln_type = "explosive"
    end

    # Return the solution

    soln = FirstOrderSolutionDet(steady_state[1:nx],hx,steady_state[nx+1:nv],gx,grc,soln_type)
    return soln

end

"""
Implements Klein's (2000) method to solve a stochastic model to first-order accuracy.

Internal function; not exposed to users.
"""
function solve_first_order_stoch(model::REModel,scheme::PerturbationScheme)

    nx = model.number_states
    ny = model.number_jumps
    nv = model.number_variables
    ns = model.number_shocks

    steady_state = scheme.steady_state
    cutoff = scheme.cutoff

    T = eltype(steady_state)

    # Compute the first-order deterministic terms

    d = compute_linearization(model,steady_state)

    @views a = d[1:nv,1:nv]
    @views b = -d[1:nv,nv+1:2*nv]

    r = schur(complex(a),complex(b))

    # Reorder the generalized eigenvalues so that those with modulus greater
    # than "cutoff" reside at the bottom.

    sel = (abs.(diag(r.S)./diag(r.T)) .< cutoff)
    ordschur!(r,sel)

    s = r.S
    t = r.T
    #q = r.Q'  # So now q*a*z = s and q*b*z = t
    z = r.Z

    # Construct the rational expectations equilibrium by eliminating the
    # unstable dynamics

    @views z11 = z[1:nx,1:nx]
    @views z21 = z[nx+1:nv,1:nx]
    @views t11 = t[1:nx,1:nx]
    @views s11 = s[1:nx,1:nx]

    if rank(z11) != nx
        error("The model does not have a stable solution.")
    end

    hx = real((z11/t11)*(s11/z11))
    gx = real(z21/z11)

    # Compute the first-order stochastic terms

    @views k = -d[1:nx,2*nv+1:2*nv+ns]
    sigma = eye(T,ns)

    # Calculate the number of unstable eigenvalues

    grc = nv - sum(sel)

    soln_type = "determinate"
    if grc < ny
        soln_type = "indeterminate"
    elseif grc > ny
        soln_type = "explosive"
    end

    # Return the solution

    soln = FirstOrderSolutionStoch(steady_state[1:nx],hx,k,steady_state[nx+1:nv],gx,sigma,grc,soln_type)
    return soln

end

"""
Implements Gomme and Klein's (2011) method to solve a model to second-order accuracy.

Exported function.
"""
function solve_second_order(model::REModel,scheme::PerturbationScheme)

    if scheme.order != "second"
        error("A second order perturbation must be specified.")
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

"""
Implements Gomme and Klein's (2011) method to solve a deterministic model to second-order accuracy.

Internal function; not exposed to users.
"""
function solve_second_order_det(model::REModel,scheme::PerturbationScheme)

    nx = model.number_states
    ny = model.number_jumps
    nv = model.number_variables
    ne = model.number_equations

    steady_state = scheme.steady_state
    cutoff = scheme.cutoff

    T = eltype(steady_state)

    # Calculate the first derivatives at the steady state

    d = compute_linearization(model,steady_state)

    # Calculate the Hessian at the steady state for each equation and stack the
    # Hessians vertically.

    point = [steady_state; steady_state]
    deriv2 = zeros(ne*2*nv,2*nv)
    for i = 1:ne
        deriv2[(i-1)*2*nv+1:i*2*nv,:] = ForwardDiff.hessian(model.each_eqn_function[i],point,ForwardDiff.HessianConfig(model.each_eqn_function[i],point,ForwardDiff.Chunk{2}()))[1:2*nv,1:2*nv]
    end

    # Compute the first-order solution

    first_order_soln = solve_first_order(model,PerturbationScheme(steady_state,cutoff,"first"))
    hx = first_order_soln.hx
    gx = first_order_soln.gx
    grc = first_order_soln.grc
    soln_type = first_order_soln.soln_type

    # Compute the second-order deterministic terms

    # Construct partitioned first derivative matrices

    @views fx  = d[:,1:nx]
    @views fy  = d[:,(nx+1):nv]
    @views fxp = d[:,(nv+1):(nv+nx)]
    @views fyp = d[:,(nv+nx+1):2*nv]

    # Set up the Sylvester equation needed to construct the matrices on the
    # second-order states.

    Mx = [I; gx; hx; gx*hx]
    q  = kron_prod_times_matrix(eye(T,nv),Mx',deriv2*Mx) #(eye(T,nv) ⊗ m')*(deriv2*m)
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
    hxx = Matrix(reshape(z[1:nx*nx,:]',nx*nx,nx)')
    gxx = Matrix(reshape(z[nx^2+1:end,:]',nx*nx,ny)')

    # Return the solution

    soln = SecondOrderSolutionDet(steady_state[1:nx],hx,hxx,steady_state[nx+1:nv],gx,gxx,grc,soln_type)
    return soln

end

"""
Implements Gomme and Klein's (2011) method to solve a stochastic model to second-order accuracy.

Internal function; not exposed to users.
"""
function solve_second_order_stoch(model::REModel, scheme::PerturbationScheme)

    nx = model.number_states
    ny = model.number_jumps
    nv = model.number_variables
    ne = model.number_equations
    ns = model.number_shocks

    steady_state = scheme.steady_state
    cutoff = scheme.cutoff

    T = eltype(steady_state)

    # Calculate the first derivatives at the steady state

    d = compute_linearization(model,steady_state)

    # Calculate the Hessian at the steady state for each equation and stack the
    # Hessians vertically.

    point = [steady_state; steady_state; zeros(ns)]
    deriv2 = zeros(ne*2*nv,2*nv)
    for i = 1:ne
        deriv2[(i-1)*2*nv+1:i*2*nv,:] = ForwardDiff.hessian(model.each_eqn_function[i], point, ForwardDiff.HessianConfig(model.each_eqn_function[i], point, ForwardDiff.Chunk{2}()))[1:2*nv,1:2*nv]
    end

    # Compute the first-order solution 

    first_order_soln = solve_first_order(model, PerturbationScheme(steady_state, cutoff, "first"))
    hx        = first_order_soln.hx
    gx        = first_order_soln.gx
    grc       = first_order_soln.grc
    soln_type = first_order_soln.soln_type
    k         = first_order_soln.k
    sigma     = first_order_soln.sigma

    # Compute the second-order deterministic terms

    # Construct partitioned first derivative matrices

    @views fx  = d[:,1:nx]
    @views fy  = d[:,(nx+1):nv]
    @views fxp = d[:,(nv+1):(nv+nx)]
    @views fyp = d[:,(nv+nx+1):2*nv]

    # Set up the Sylvester equation needed to construct the matrices on the
    # second-order states.

    Mx = [I; gx; hx; gx * hx]
    q = kron_prod_times_matrix(eye(T,nv),Mx',deriv2*Mx) #(eye(T,nv) ⊗ m')*(deriv2*m)
    b1 = kron(fxp,eye(T,nx)) #fxp ⊗ eye(T,nx) 
    b2 = kron(fyp,eye(T,nx)) #fyp ⊗ eye(T,nx)
    b4 = kron(fy,eye(T,nx)) #fy ⊗ eye(T,nx) 
    c1 = kron(eye(T,ny),hx') #eye(T,ny) ⊗ hx'
    c2 = kron(gx,eye(T,nx)) #gx ⊗ eye(T,nx)

    # Use a Sylvester equation solver to compute the second order terms on the states.

    A = [b1+b2*c2 b4]
    B = [zeros(size(b2,1),nx*nx) b2*c1]
    C = hx
    D = q

    B .= A\B
    D .= A\D

    z = dsylvester(B,C,-D)
    hxx = z[1:nx*nx,:]
    gxx = z[nx^2+1:end,:]

    # Compute the second-order stochastic terms

    # Set up the LP problem needed to construct the intercepts, which contain
    # the volatility effects

    qq = [fxp+fyp*gx fyp+fy]
    nn = [zeros(nv,nx); I; gx]
    q = fyp*tracem(kron_prod_times_matrix(eye(T,ny),k*sigma'*k',gxx)) + tracem(kron_prod_times_matrix(eye(T,nv),nn',deriv2*nn*k*sigma'*k'))
    ss = -qq\q
    hss = ss[1:nx]
    gss = ss[nx+1:nv]

    # Reshape hxx and gxx to be conformable with Binning, Levintal, etc.

    hxx = Matrix(reshape(z[1:nx*nx,:]',nx*nx,nx)')
    gxx = Matrix(reshape(z[nx^2+1:end,:]',nx*nx,ny)')

    # Return the solution

    soln = SecondOrderSolutionStoch(steady_state[1:nx],hx,hss,hxx,k,steady_state[nx+1:nv],gx,gss,gxx,sigma,grc,soln_type)
    return soln

end

"""
Implements Binning's (2013) method to solve a model to third-order accuracy.

Exported function.
"""
function solve_third_order(model::REModel,scheme::PerturbationScheme,skewness::Union{Array{T,1},Array{T,2}} = zeros(model.number_shocks,model.number_shocks^2)) where {T<:Real}

    if scheme.order != "third"
        error("A third order perturbation must be supplied.")
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

"""
Implements Binning's (2013) method to solve a deterministic model to third-order accuracy.

Internal function; not exposed to users.
"""
function solve_third_order_det(model::REModel,scheme::PerturbationScheme)

    nv = model.number_variables
    nx = model.number_states
    ny = model.number_jumps
    ne = model.number_equations

    steady_state = scheme.steady_state
    cutoff = scheme.cutoff

    # Compute the first, second, and third derivatives

    model_equations = model.each_eqn_function
    point = [steady_state; steady_state]

    first_derivs  = zeros(ne,2*nv)
    second_derivs = zeros(ne,4*nv^2)
    third_derivs  = zeros(ne,8*nv^3)

    for i = 1:ne

        first_d(x) = ForwardDiff.gradient(model_equations[i],x,ForwardDiff.GradientConfig(model_equations[i],x,ForwardDiff.Chunk{2}()))[1:2*nv]
        first_derivs[i,:] .= first_d(point)

        #second_d(x) = ForwardDiff.hessian(model_equations[i],x,ForwardDiff.HessianConfig(model_equations[i],x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
        second_d(x) = ForwardDiff.jacobian(first_d,x,ForwardDiff.JacobianConfig(first_d,x,ForwardDiff.Chunk{1}()))[:,1:2*nv]
        second_derivs[i,:] .= vec(second_d(point))

        third_d(x) = ForwardDiff.jacobian(second_d,x,ForwardDiff.JacobianConfig(second_d,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
        third_derivs[i,:] .= vec(third_d(point))

    end

    # Compute the first-order solution

    first_order_soln = solve_first_order(model,PerturbationScheme(steady_state,cutoff,"first"))
    hx        = first_order_soln.hx
    gx        = first_order_soln.gx
    grc       = first_order_soln.grc
    soln_type = first_order_soln.soln_type

    # Compute the second-order deterministic terms

    Mx = [I; gx; hx; gx*hx] # Gradient of the policy functions

    @views fx  = first_derivs[:,1:nx]
    @views fy  = first_derivs[:,nx+1:nv]
    @views fxp = first_derivs[:,nv+1:nv+nx]
    @views fyp = first_derivs[:,nv+nx+1:2*nv]

    A = [fxp+fyp*gx fy]
    B = [zeros(nv,nx) fyp]
    D = -matrix_times_kron_prod(second_derivs,Mx,Mx)

    z = martin_van_loan(A,B,hx,D,1)
    hxx = z[1:nx,:]
    gxx = z[nx+1:nv,:]

    # Compute the third-order deterministic terms

    Mxx = [zeros(nx,nx^2); gxx; hxx; gxx*kron(hx,hx) + gx*hxx] # Hessian of the policy functions

    omega1 = create_omega3(nx)

    D = -matrix_times_kron_prod(third_derivs,Mx,Mx,Mx) - matrix_times_kron_prod(second_derivs,Mx,Mxx)*omega1 - fyp*gxx*kron(hx,hxx)*omega1

    z = martin_van_loan(A,B,hx,D,2)
    hxxx = z[1:nx,:]
    gxxx = z[nx+1:nv,:]

    # Return the solution

    soln = ThirdOrderSolutionDet(steady_state[1:nx],hx,hxx,hxxx,steady_state[nx+1:nv],gx,gxx,gxxx,grc,soln_type)
    return soln

end

"""
Implements Binning's (2013) method to solve a stochastic model to third-order accuracy.

Internal function; not exposed to users.
"""
function solve_third_order_stoch(model::REModel,scheme::PerturbationScheme,skewness::Union{Array{T,1},Array{T,2}}) where {T<:Real}

    ns = model.number_shocks
    nv = model.number_variables
    nx = model.number_states
    ny = model.number_jumps
    ne = model.number_equations

    steady_state = scheme.steady_state
    cutoff = scheme.cutoff

    # Compute the first, second, and third derivatives

    model_equations = model.each_eqn_function
    point = [steady_state; steady_state; zeros(ns)]

    first_derivs  = zeros(ne,2*nv)
    second_derivs = zeros(ne,4*nv^2)
    third_derivs  = zeros(ne,8*nv^3)

    for i = 1:ne

        first_d(x) = ForwardDiff.gradient(model_equations[i],x,ForwardDiff.GradientConfig(model_equations[i],x,ForwardDiff.Chunk{2}()))[1:2*nv]
        first_derivs[i,:] .= first_d(point)

        #second_d(x) = ForwardDiff.hessian(model_equations[i],x,ForwardDiff.HessianConfig(model_equations[i],x,ForwardDiff.Chunk{2}()))[1:2*nv,1:2*nv]
        second_d(x) = ForwardDiff.jacobian(first_d,x,ForwardDiff.JacobianConfig(first_d,x,ForwardDiff.Chunk{1}()))[1:2*nv,1:2*nv]
        second_derivs[i,:] .= vec(second_d(point))

        third_d(x) = ForwardDiff.jacobian(second_d,x,ForwardDiff.JacobianConfig(second_d,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
        third_derivs[i,:] .= vec(third_d(point))

    end

    # Compute the first-order solution

    first_order_soln = solve_first_order(model,PerturbationScheme(steady_state,cutoff,"first"))
    hx        = first_order_soln.hx
    gx        = first_order_soln.gx
    grc       = first_order_soln.grc
    soln_type = first_order_soln.soln_type
    k         = first_order_soln.k
    sigma     = first_order_soln.sigma

    # Compute the second-order deterministic terms

    Mx = [I; gx; hx; gx*hx] # Gradient of the policy functions

    @views fx  = first_derivs[:,1:nx]
    @views fy  = first_derivs[:,nx+1:nv]
    @views fxp = first_derivs[:,nv+1:nv+nx]
    @views fyp = first_derivs[:,nv+nx+1:2*nv]

    A = [fxp+fyp*gx fy]
    B = [zeros(nv,nx) fyp]
    D = -matrix_times_kron_prod(second_derivs,Mx,Mx)

    z = martin_van_loan(A,B,hx,D,1)
    hxx = z[1:nx,:]
    gxx = z[nx+1:nv,:]

    # Compute the second-order stochastic terms

    Ns = [zeros(nv,nx); I; gx]

    G = [fxp + fyp*gx fy + fyp]
    F = -trm(matrix_times_kron_prod(second_derivs,Ns,Ns*k*sigma*k')) .- fyp*trm(matrix_times_kron_prod(gxx,eye(T,nx),k*sigma*k'))

    z = G\F
    hss = z[1:nx]
    gss = z[nx+1:nv]

    # Compute the third-order deterministic terms

    Mxx = [zeros(nx,nx^2); gxx; hxx; gxx*kron(hx,hx) + gx*hxx] # Hessian of the policy functions

    omega1 = create_omega3(nx)

    D = -matrix_times_kron_prod(third_derivs,Mx,Mx,Mx) .- matrix_times_kron_prod(second_derivs,Mx,Mxx)*omega1 .- fyp*gxx*kron(hx,hxx)*omega1

    z = martin_van_loan(A,B,hx,D,2)
    hxxx = z[1:nx,:]
    gxxx = z[nx+1:nv,:]

    # Compute the third-order stochastic terms

    Msx = [zeros(nv + nx,nx^2); gxx*kron(hx,eye(T,nx))]
    Mss = [zeros(nx,1); gss; hss; gx*hss + trm(gxx*kron(eye(T,nx),k*sigma*k')) + gss]

    C = hx
    F = (-trm2(matrix_times_kron_prod(third_derivs,Mx,Ns,Ns*k*sigma*k')) .- 2*trm2(matrix_times_kron_prod(second_derivs,Msx,Ns*k*sigma*k')) .- matrix_times_kron_prod(second_derivs,Mx,Mss[:,:])
        .-fyp*(trm2(gxxx*kron(hx,eye(T,nx^2))*kron(eye(T,nx^2),k*sigma*k')) .+ matrix_times_kron_prod(gxx,hx,hss[:,:])))

    F .= A\F

    z = dsylvester(B,C,F)
    hssx = z[1:nx,:]
    gssx = z[nx+1:nv,:]

    # Compute the third-order stochastic terms related to skewness

    hsss = zeros(nx)
    gsss = zeros(ny)

    Nss = [zeros(nv+nx,nx^2); gxx]

    if sum(isequal.(skewness,0.0)) == 0

        skew = zeros(nx,nx^2)
        skew[1:ns,1:ns^2] = skewness

        F = -trm(third_derivs*kron(Ns,kron(Ns,Ns*skew))) .- 3*trm(second_derivs*kron(Nss,Ns*skew)) .- fyp*trm(gxxx*kron(eye(T,nx^2),skew))
        z = G\F
        hsss .= z[1:nx]
        gsss .= z[nx+1:nv]

    end

    # Return the solution

    soln = ThirdOrderSolutionStoch(steady_state[1:nx],hx,hss,hxx,hsss,hssx,hxxx,k,steady_state[nx+1:nv],gx,gss,gxx,gsss,gssx,gxxx,sigma,skewness,grc,soln_type)
    return soln

end

"""
Solves models to fourth-order accuracy, drawing on Binning (2013) and Levintal (2017).

Exported function.
"""
function solve_fourth_order(model::REModel,scheme::PerturbationScheme)

    if scheme.order != "fourth"
        error("A fourth order perturbation must be supplied.")
    end

    ns = model.number_shocks
    if ns == 0
        soln = solve_fourth_order_det(model,scheme)
        return soln
    else
        soln = solve_fourth_order_stoch(model,scheme)
        return soln
    end

end

"""
Solves deterministic models to fourth-order accuracy, drawing on Binning (2013) and Levintal (2017).

Internal function; not exposed to users.
"""
function solve_fourth_order_det(model::REModel,scheme::PerturbationScheme)

    nv = model.number_variables
    nx = model.number_states
    ny = model.number_jumps
    ne = model.number_equations

    steady_state = scheme.steady_state
    cutoff = scheme.cutoff

    # Compute the first, second, third, and fourth derivatives

    model_equations = model.each_eqn_function
    point = [steady_state; steady_state]

    first_derivs  = zeros(ne,2*nv)
    second_derivs = zeros(ne,4*nv^2)
    third_derivs  = zeros(ne,8*nv^3)
    fourth_derivs = zeros(ne,16*nv^4)

    for i = 1:ne

        first_d(x) = ForwardDiff.gradient(model_equations[i],x,ForwardDiff.GradientConfig(model_equations[i],x,ForwardDiff.Chunk{1}()))[1:2*nv]
        first_derivs[i,:] .= first_d(point)

        #second_d(x) = ForwardDiff.hessian(model_equations[i],x,ForwardDiff.HessianConfig(model_equations[i],x,ForwardDiff.Chunk{1}()))[:,1:2*nv]
        second_d(x) = ForwardDiff.jacobian(first_d,x,ForwardDiff.JacobianConfig(first_d,x,ForwardDiff.Chunk{1}()))[1:2*nv,1:2*nv]
        second_derivs[i,:] .= vec(second_d(point))

        third_d(x) = ForwardDiff.jacobian(second_d,x,ForwardDiff.JacobianConfig(second_d,x,ForwardDiff.Chunk{1}()))[:,1:2*nv]
        third_derivs[i,:] .= vec(third_d(point))

        fourth_d(x) = ForwardDiff.jacobian(third_d,x,ForwardDiff.JacobianConfig(third_d,x,ForwardDiff.Chunk{1}()))[:,1:2*nv]
        fourth_derivs[i,:] .= vec(fourth_d(point))

    end

    # Compute the first-order solution

    first_order_soln = solve_first_order(model,PerturbationScheme(steady_state,cutoff,"first"))
    hx        = first_order_soln.hx
    gx        = first_order_soln.gx
    grc       = first_order_soln.grc
    soln_type = first_order_soln.soln_type

    # Compute the second-order deterministic terms

    Mx = [I; gx; hx; gx*hx] # Gradient of the policy functions

    @views fx  = first_derivs[:,1:nx]
    @views fy  = first_derivs[:,nx+1:nv]
    @views fxp = first_derivs[:,nv+1:nv+nx]
    @views fyp = first_derivs[:,nv+nx+1:2*nv]

    A = [fxp+fyp*gx fy]
    B = [zeros(nv,nx) fyp]
    D = -matrix_times_kron_prod(second_derivs,Mx,Mx)

    z = martin_van_loan(A,B,hx,D,1)
    hxx = z[1:nx,:]
    gxx = z[nx+1:nv,:]

    # Compute the third-order deterministic terms

    Mxx = [zeros(nx,nx^2); gxx; hxx; gxx*kron(hx,hx) + gx*hxx] # Hessian of the policy functions

    omega1 = create_omega3(nx)

    D = -matrix_times_kron_prod(third_derivs,Mx,Mx,Mx) - matrix_times_kron_prod(second_derivs,Mx,Mxx)*omega1 - fyp*gxx*kron(hx,hxx)*omega1

    z = martin_van_loan(A,B,hx,D,2)
    hxxx = z[1:nx,:]
    gxxx = z[nx+1:nv,:]

    # Compute the fourth-order deterministic terms

    Mxxx = [zeros(nx,nx^3); gxxx; hxxx; gxxx*kron(hx,hx,hx) + gxx*kron(hxx,hx)*omega1 + gx*hxxx] # Third derivatives of the policy functions

    omega2, omega3, omega4 = create_omega4(nx)

    D = (-matrix_times_kron_prod(fourth_derivs,[Mx, Mx, Mx, Mx]) - matrix_times_kron_prod(third_derivs,[Mx, Mx, Mxx])*omega2 - matrix_times_kron_prod(second_derivs,Mx,Mxxx)*omega3 - matrix_times_kron_prod(second_derivs,Mxx,Mxx)*omega4
        -fyp*(matrix_times_kron_prod(gxxx,[hx, hx, hxx])*omega2 + matrix_times_kron_prod(gxx,hx,hxxx)*omega3 + matrix_times_kron_prod(gxx,hxx,hxx)*omega4))

    z = martin_van_loan(A,B,hx,D,3)
    hxxxx = z[1:nx,:]
    gxxxx = z[nx+1:nv,:]

    # Return the solution

    soln = FourthOrderSolutionDet(steady_state[1:nx],hx,hxx,hxxx,hxxxx,steady_state[nx+1:nv],gx,gxx,gxxx,gxxxx,grc,soln_type)
    return soln

end

"""
Solves stochastic models to fourth-order accuracy, drawing on Binning (2013) and Levintal (2017).

Internal function; not exposed to users.
"""
function solve_fourth_order_stoch(model::REModel,scheme::PerturbationScheme)

    ns = model.number_shocks
    nv = model.number_variables
    nx = model.number_states
    ny = model.number_jumps
    ne = model.number_equations

    steady_state = scheme.steady_state
    cutoff = scheme.cutoff

    # Compute the first, second, third and fourth derivatives

    model_equations = model.each_eqn_function
    point = [steady_state; steady_state; zeros(ns)]

    first_derivs  = zeros(ne,2*nv)
    second_derivs = zeros(ne,4*nv^2)
    third_derivs  = zeros(ne,8*nv^3)
    fourth_derivs = zeros(ne,16*nv^4)

    for i = 1:ne

        first_d(x) = ForwardDiff.gradient(model_equations[i],x,ForwardDiff.GradientConfig(model_equations[i],x,ForwardDiff.Chunk{1}()))[1:2*nv]
        first_derivs[i,:] .= first_d(point)

        #second_d(x) = ForwardDiff.hessian(model_equations[i],x,ForwardDiff.HessianConfig(model_equations[i],x,ForwardDiff.Chunk{1}()))[1:2*nv,1:2*nv]
        second_d(x) = ForwardDiff.jacobian(first_d,x,ForwardDiff.JacobianConfig(first_d,x,ForwardDiff.Chunk{1}()))[1:2*nv,1:2*nv]
        second_derivs[i,:] .= vec(second_d(point))

        third_d(x) = ForwardDiff.jacobian(second_d,x,ForwardDiff.JacobianConfig(second_d,x,ForwardDiff.Chunk{1}()))[:,1:2*nv]
        third_derivs[i,:] .= vec(third_d(point))

        fourth_d(x) = ForwardDiff.jacobian(third_d,x,ForwardDiff.JacobianConfig(third_d,x,ForwardDiff.Chunk{1}()))[:,1:2*nv]
        fourth_derivs[i,:] .= vec(fourth_d(point))

    end

    # Compute the first-order Solution

    first_order_soln = solve_first_order(model,PerturbationScheme(steady_state,cutoff,"first"))
    hx        = first_order_soln.hx
    gx        = first_order_soln.gx
    grc       = first_order_soln.grc
    soln_type = first_order_soln.soln_type
    k         = first_order_soln.k
    sigma     = first_order_soln.sigma

    # Compute the second-order deterministic terms

    Mx = [I; gx; hx; gx*hx] # Gradient of the policy functions

    @views fx  = first_derivs[:,1:nx]
    @views fy  = first_derivs[:,nx+1:nv]
    @views fxp = first_derivs[:,nv+1:nv+nx]
    @views fyp = first_derivs[:,nv+nx+1:2*nv]

    A = [fxp+fyp*gx fy]
    B = [zeros(nv,nx) fyp]
    D = -matrix_times_kron_prod(second_derivs,Mx,Mx)

    z = martin_van_loan(A,B,hx,D,1)
    hxx = z[1:nx,:]
    gxx = z[nx+1:nv,:]

    # Compute the second-order stochastic terms

    Ns = [zeros(nv,nx); I; gx]

    G = [fxp+fyp*gx fy+fyp]
    F = -trm(matrix_times_kron_prod(second_derivs,Ns,Ns*k*sigma*k')) .- fyp*trm(matrix_times_kron_prod(gxx,eye(nx),k*sigma*k'))

    z = G\F
    hss = z[1:nx]
    gss = z[nx+1:nv]

    # Compute the third-order deterministic terms

    Mxx = [zeros(nx,nx^2); gxx; hxx; gxx*kron(hx,hx) + gx*hxx] # Hessian of the policy functions

    omega1 = create_omega3(nx)

    D = -matrix_times_kron_prod(third_derivs,Mx,Mx,Mx) .- matrix_times_kron_prod(second_derivs,Mx,Mxx)*omega1 .- fyp*gxx*kron(hx,hxx)*omega1

    z = martin_van_loan(A,B,hx,D,2)
    hxxx = z[1:nx,:]
    gxxx = z[nx+1:nv,:]

    # Compute the third-order stochastic terms

    Mxs = [zeros(nv+nx,nx^2); gxx*kron(hx,eye(nx))]
    Mss = [zeros(nx,1); gss; hss; gx*hss+trm(gxx*kron(eye(nx),k*sigma*k'))+gss]

    C = hx
    F = (-trm2(matrix_times_kron_prod(third_derivs,Mx,Ns,Ns*k*sigma*k')) .- 2*trm2(matrix_times_kron_prod(second_derivs,Mxs,Ns*k*sigma*k')) .- matrix_times_kron_prod(second_derivs,Mx,Mss[:,:])
        .-fyp*(trm2(gxxx*kron(hx,eye(nx^2))*kron(eye(nx^2),k*sigma*k')) .+ matrix_times_kron_prod(gxx,hx,hss[:,:])))

    F .= A\F

    z = dsylvester(B,C,F)
    hssx = z[1:nx,:]
    gssx = z[nx+1:nv,:]

    #hsxx = zeros(nx,nx^2)
    #gsxx = zeros(ny,nx^2)

    #hsss = zeros(nx)
    #gsss = zeros(ny)

    # Compute the fourth-order deterministic terms

    Mxxx = [zeros(nx,nx^3); gxxx; hxxx; gxxx*kron(hx,hx,hx) + gxx*kron(hx,hxx)*omega1 + gx*hxxx] # Third derivatives of the policy functions

    omega2, omega3, omega4 = create_omega4(nx)

    D = (-matrix_times_kron_prod(fourth_derivs,[Mx,Mx,Mx,Mx]) .- matrix_times_kron_prod(third_derivs,[Mx,Mx,Mxx])*omega2 .- matrix_times_kron_prod(second_derivs,Mx,Mxxx)*omega3 .- matrix_times_kron_prod(second_derivs,Mxx,Mxx)*omega4
        .-fyp*(matrix_times_kron_prod(gxxx,[hx,hx,hxx])*omega2 .+ matrix_times_kron_prod(gxx,hx,hxxx)*omega3 .+ matrix_times_kron_prod(gxx,hxx,hxx)*omega4))

    z = martin_van_loan(A,B,hx,D,3)
    hxxxx = z[1:nx,:]
    gxxxx = z[nx+1:nv,:]

    # Compute the fourth-order stochastic terms

    Mxxs = zeros(2*nx+2*ny,nx^2)
    Mxss = [zeros(nx,nx); gssx; hssx; gx*hssx + gxx*(kron(hss,hx)) + trm2(matrix_times_kron_prod(gxxx,eye(nx),k*sigma*k',hx))]
    Msss = zeros(2*nx+2*ny,1)

    D = (- matrix_times_kron_prod(second_derivs,Mss[:,:],Mxx) .- matrix_times_kron_prod(second_derivs,Mxss,Mx) .- matrix_times_kron_prod(second_derivs,Mx,Mxss) .- matrix_times_kron_prod(third_derivs,[Mss[:,:],Mx,Mx]) .- trm3(matrix_times_kron_prod(fourth_derivs,[Ns,Ns*k*sigma*k',Mx,Mx])) .- trm3(matrix_times_kron_prod(third_derivs,[Ns,Ns*k*sigma*k',Mxx]))
        .- fyp*(matrix_times_kron_prod(gxx,hss[:,:],hxx) .+ gssx*hxx + gxx*kron(hssx,hx) + gxx*kron(hx,hssx) .+ matrix_times_kron_prod(gxxx,[hss[:,:],hx,hx]) .+ trm3(matrix_times_kron_prod(gxxxx,[eye(nx),k*sigma*k',hx,hx]))))
    
    z = martin_van_loan(A,B,hx,D,1)
    hxxss = z[1:nx,:]
    gxxss = z[nx+1:nv,:]

    G = [fxp+fyp*gx fy+fyp]

    F = (- 3*matrix_times_kron_prod(second_derivs,Mss,Mss) .- fyp*(3*matrix_times_kron_prod(gxx,hss[:,:],hss[:,:]) .+ 3*trm(matrix_times_kron_prod(gxxx,[eye(nx),k*sigma*k',hss[:,:]])) .+ gssx*hss))

    z = G\F
    hssss = z[1:nx]
    gssss = z[nx+1:nv]

    # Return the solution
    
    soln = FourthOrderSolutionStoch(steady_state[1:nx],hx,hss,hxx,hssx,hxxx,hssss,hxxss,hxxxx,k,steady_state[nx+1:nv],gx,gss,gxx,gssx,gxxx,gssss,gxxss,gxxxx,sigma,grc,soln_type)
    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{ChebyshevSchemeDet,ChebyshevSchemeOBCDet})

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number = scheme.node_number
    order = scheme.order
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: ChebyshevSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.ftol)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],Tuple(length.(grid)))
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],Tuple(length.(grid)))
    end

    if typeof(order) <: Integer
        ord = Tuple(fill(order,nx))
    else
        ord = Tuple(order)
    end

    weights = [zeros(ord.+1) for _ in 1:length(jumps_approximated)]

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    state = Array{T,1}(undef,nx)
    init = Array{T,1}(undef,nv)

    N = prod(length.(grid))

    if node_generator == chebyshev_nodes
        cheb_weights = chebyshev_weights
    elseif node_generator == chebyshev_extrema
        cheb_weights = chebyshev_weights_extrema
    elseif node_generator == chebyshev_extended
        cheb_weights = chebyshev_weights_extended
    end

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= cheb_weights(variables[jumps_approximated[i]],grid,order,domain)
        end

        for i = 1:N

            sub = ind2sub(i,Tuple(length.(grid)))
            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_chebyshev(state,weights,order,domain)
            if typeof(scheme) <: ChebyshevSchemeOBCDet
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,maxiters, iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionDet([variables[ny+1:end]; variables[1:ny]],weights,grid,order,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{ChebyshevSchemeDet,ChebyshevSchemeOBCDet},threads::S) where {S<:Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number = scheme.node_number
    order = scheme.order
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: ChebyshevSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.ftol)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],Tuple(length.(grid)))
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],Tuple(length.(grid)))
    end

    if typeof(order) <: Integer
        ord = Tuple(fill(order,nx))
    else
        ord = Tuple(order)
    end

    weights = [zeros(ord.+1) for _ in 1:length(jumps_approximated)]

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    N = prod(length.(grid))

    if node_generator == chebyshev_nodes
        cheb_weights = chebyshev_weights
    elseif node_generator == chebyshev_extrema
        cheb_weights = chebyshev_weights_extrema
    elseif node_generator == chebyshev_extended
        cheb_weights = chebyshev_weights_extended
    end

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= cheb_weights(variables[jumps_approximated[i]],grid,order,domain)
        end

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N

                sub = ind2sub(i,Tuple(length.(grid)))
                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_chebyshev(state,weights,order,domain)
                if typeof(scheme) <: ChebyshevSchemeOBCDet
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionDet([variables[ny+1:end]; variables[1:ny]],weights,grid,order,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{ChebyshevSchemeStoch,ChebyshevSchemeOBCStoch})

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    order = scheme.order
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: ChebyshevSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.ftol)

    d = compute_linearization(model,initial_guess)
    GAMMA = d[1:ns,nv+1:nv+ns]
    k   = -GAMMA\d[1:ns,2*nv+1:end]
    RHO = -GAMMA\d[1:ns,1:ns]
    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    (eps_nodes, eps_weights) = hermite(num_quad_nodes)
    integrals = compute_chebyshev_integrals(eps_nodes,eps_weights,grid,order,RHO,k)

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],Tuple(length.(grid)))
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],Tuple(length.(grid)))
    end

    if typeof(order) <: Integer
        ord = Tuple(fill(order,nx))
    else
        ord = Tuple(order)
    end

    weights = [zeros(ord.+1) for _ in 1:length(jumps_approximated)]
    scaled_weights = [zeros(ord.+1) for _ in 1:length(jumps_approximated)]

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    state = Array{T,1}(undef,nx)
    init = Array{T,1}(undef,nv)

    N = prod(length.(grid))

    if node_generator == chebyshev_nodes
        cheb_weights = chebyshev_weights
    elseif node_generator == chebyshev_extrema
        cheb_weights = chebyshev_weights_extrema
    elseif node_generator == chebyshev_extended
        cheb_weights = chebyshev_weights_extended
    end

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= cheb_weights(variables[jumps_approximated[i]],grid,order,domain)
        end

        scale_chebyshev_weights!(weights,scaled_weights,integrals,jumps_approximated,ns)

        for i = 1:N

            sub = ind2sub(i,Tuple(length.(grid)))
            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_chebyshev(state,scaled_weights,order,domain)
            if typeof(scheme) <: ChebyshevSchemeOBCStoch
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionStoch([variables[ny+1:end]; variables[1:ny]],weights,integrals,grid,order,domain,k,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{ChebyshevSchemeStoch,ChebyshevSchemeOBCStoch},threads::S) where {S<:Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    order = scheme.order
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: ChebyshevSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.ftol)

    d = compute_linearization(model,initial_guess)
    GAMMA = d[1:ns,nv+1:nv+ns]
    k   = -GAMMA\d[1:ns,2*nv+1:end]
    RHO = -GAMMA\d[1:ns,1:ns]
    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = compute_chebyshev_integrals(eps_nodes,eps_weights,grid,order,RHO,k)

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],Tuple(length.(grid)))
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],Tuple(length.(grid)))
    end

    if typeof(order) <: Integer
        ord = Tuple(fill(order,nx))
    else
        ord = Tuple(order)
    end

    weights = [zeros(ord.+1) for _ in 1:length(jumps_approximated)]
    scaled_weights = [zeros(ord.+1) for _ in 1:length(jumps_approximated)]

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    N = prod(length.(grid))

    if node_generator == chebyshev_nodes
        cheb_weights = chebyshev_weights
    elseif node_generator == chebyshev_extrema
        cheb_weights = chebyshev_weights_extrema
    elseif node_generator == chebyshev_extended
        cheb_weights = chebyshev_weights_extended
    end

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= cheb_weights(variables[jumps_approximated[i]],grid,order,domain)
        end

        scale_chebyshev_weights!(weights,scaled_weights,integrals,jumps_approximated,ns)

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N

                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end

                init = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_chebyshev(state,scaled_weights,order,domain)
                if typeof(scheme) <: ChebyshevSchemeOBCStoch
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end
    
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionStoch([variables[ny+1:end]; variables[1:ny]],weights,integrals,grid,order,domain,k,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{ChebyshevSchemeDet,ChebyshevSchemeOBCDet}) where {R<:Union{PerturbationSolutionDet,ProjectionSolutionDet}}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    node_generator = scheme.node_generator
    node_number = scheme.node_number
    order = scheme.order
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: ChebyshevSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    ss_eqm = state_space_eqm(soln)

    T = typeof(scheme.ftol)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    N = prod(length.(grid))
    state = Array{T,1}(undef,nx)

    for j = 1:N
        sub = ind2sub(j,Tuple(length.(grid)))

        for i = 1:nx
            state[i] = grid[i][sub[i]]
        end

        dr = ss_eqm.g(state)
        te = ss_eqm.h(state)

        for i = 1:ny
            variables[i][j] = min(max(dr[i],lb[i]),ub[i])
        end
        for i = 1:nx
            variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
        end

    end

    if typeof(order) <: Integer
        ord = Tuple(fill(order,nx))
    else
        ord = Tuple(order)
    end

    weights = [zeros(ord.+1) for _ in 1:length(jumps_approximated)]

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    init = Array{T,1}(undef,nv)

    if node_generator == chebyshev_nodes
        cheb_weights = chebyshev_weights
    elseif node_generator == chebyshev_extrema
        cheb_weights = chebyshev_weights_extrema
    elseif node_generator == chebyshev_extended
        cheb_weights = chebyshev_weights_extended
    end

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= cheb_weights(variables[jumps_approximated[i]],grid,order,domain)
        end

        for i = 1:N

            sub = ind2sub(i,Tuple(length.(grid)))
            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_chebyshev(state,weights,order,domain)
            if typeof(scheme) <: ChebyshevSchemeOBCDet
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionDet([variables[ny+1:end]; variables[1:ny]],weights,grid,order,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{ChebyshevSchemeDet,ChebyshevSchemeOBCDet},threads::S) where {R<:Union{PerturbationSolutionDet,ProjectionSolutionDet},S<:Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    node_generator = scheme.node_generator
    node_number = scheme.node_number
    order = scheme.order
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: ChebyshevSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    ss_eqm = state_space_eqm(soln)

    T = typeof(scheme.ftol)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    N = prod(length.(grid))

    @sync Threads.@threads for t = 1:threads
        for j = t:threads:N
            sub = ind2sub(j,Tuple(length.(grid)))

            state = Array{T,1}(undef,nx)
            for i = 1:nx
                state[i] = grid[i][sub[i]]
            end

            dr = ss_eqm.g(state)
            te = ss_eqm.h(state)

            for i = 1:ny
                variables[i][j] = min(max(dr[i],lb[i]),ub[i])
            end
            for i = 1:nx
                variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
            end

        end
    end

    if typeof(order) <: Integer
        ord = Tuple(fill(order,nx))
    else
        ord = Tuple(order)
    end

    weights = [zeros(ord.+1) for _ in 1:length(jumps_approximated)]

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    if node_generator == chebyshev_nodes
        cheb_weights = chebyshev_weights
    elseif node_generator == chebyshev_extrema
        cheb_weights = chebyshev_weights_extrema
    elseif node_generator == chebyshev_extended
        cheb_weights = chebyshev_weights_extended
    end

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= cheb_weights(variables[jumps_approximated[i]],grid,order,domain)
        end

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N

                sub = ind2sub(i,Tuple(length.(grid)))
                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_chebyshev(state,weights,order,domain)
                if typeof(scheme) <: ChebyshevSchemeOBCDet
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionDet([variables[ny+1:end]; variables[1:ny]],weights,grid,order,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{ChebyshevSchemeStoch,ChebyshevSchemeOBCStoch}) where {R<:Union{PerturbationSolutionStoch,ProjectionSolutionStoch}}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    order = scheme.order
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: ChebyshevSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    ss_eqm = state_space_eqm(soln)

    k = soln.k[1:ns,:]

    if typeof(soln) <: PerturbationSolution
        RHO = soln.hx[1:ns,1:ns]
    elseif typeof(soln) <: ProjectionSolution
        d = compute_linearization(model,initial_guess)
        GAMMA = d[1:ns,nv+1:nv+ns]
        k   = -GAMMA\d[1:ns,2*nv+1:end]
        RHO = -GAMMA\d[1:ns,1:ns]
    end

    T = typeof(scheme.ftol)

    if domain == []
        if typeof(soln) <: ProjectionSolution
            domain = soln.domain
        elseif typeof(soln) <: PerturbationSolutionStoch
            soln_fo = FirstOrderSolutionStoch(soln.hbar,soln.hx,soln.k,soln.gbar,soln.gx,soln.sigma,soln.grc,soln.soln_type)
            state_vars,jump_vars = compute_variances(soln_fo)
            domain = Matrix([soln.hbar + 3*sqrt.(diag(state_vars)) soln.hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
        end
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    (eps_nodes, eps_weights) = hermite(num_quad_nodes)
    integrals = compute_chebyshev_integrals(eps_nodes,eps_weights,grid,order,RHO,k[1:ns,1:ns]*k[1:ns,1:ns]')

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    N = prod(length.(grid))
    state = Array{T,1}(undef,nx)

    for j = 1:N
        sub = ind2sub(j,Tuple(length.(grid)))

        for i = 1:nx
            state[i] = grid[i][sub[i]]
        end

        dr = ss_eqm.g(state)
        te = ss_eqm.h(state,zeros(ns))

        for i = 1:ny
            variables[i][j] = min(max(dr[i],lb[i]),ub[i])
        end
        for i = 1:nx
            variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
        end

    end

    if typeof(order) <: Integer
        ord = Tuple(fill(order,nx))
    else
        ord = Tuple(order)
    end

    weights = [zeros(ord.+1) for _ in 1:length(jumps_approximated)]
    scaled_weights = [zeros(ord.+1) for _ in 1:length(jumps_approximated)]

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    init = Array{T,1}(undef,nv)

    if node_generator == chebyshev_nodes
        cheb_weights = chebyshev_weights
    elseif node_generator == chebyshev_extrema
        cheb_weights = chebyshev_weights_extrema
    elseif node_generator == chebyshev_extended
        cheb_weights = chebyshev_weights_extended
    end

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= cheb_weights(variables[jumps_approximated[i]],grid,order,domain)
        end

        scale_chebyshev_weights!(weights,scaled_weights,integrals,jumps_approximated,ns)

        for i = 1:N

            sub = ind2sub(i,Tuple(length.(grid)))
            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_chebyshev(state,scaled_weights,order,domain)
            if typeof(scheme) <: ChebyshevSchemeOBCStoch
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionStoch([variables[ny+1:end]; variables[1:ny]],weights,integrals,grid,order,domain,k[1:ns,:],iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{ChebyshevSchemeStoch,ChebyshevSchemeOBCStoch},threads::S) where {R<:Union{PerturbationSolutionStoch,ProjectionSolutionStoch},S<:Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    node_number = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    order = scheme.order
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: ChebyshevSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    ss_eqm = state_space_eqm(soln)

    k = soln.k[1:ns,:]

    if typeof(soln) <: PerturbationSolution
        RHO = soln.hx[1:ns,1:ns]
    elseif typeof(soln) <: ProjectionSolution
        d = compute_linearization(model,initial_guess)
        GAMMA = d[1:ns,nv+1:nv+ns]
        k   = -GAMMA\d[1:ns,2*nv+1:end]
        RHO = -GAMMA\d[1:ns,1:ns]
    end

    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end

    T = typeof(scheme.ftol)

    if domain == []
        if typeof(soln) <: ProjectionSolution
            domain = soln.domain
        elseif typeof(soln) <: PerturbationSolutionStoch
            soln_fo = FirstOrderSolutionStoch(soln.hbar,soln.hx,soln.k,soln.gbar,soln.gx,soln.sigma,soln.grc,soln.soln_type)
            state_vars,jump_vars = compute_variances(soln_fo)
            domain = Matrix([soln.hbar + 3*sqrt.(diag(state_vars)) soln.hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
        end
    end

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = node_generator(node_number[i],domain[:,i])
    end

    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    integrals = compute_chebyshev_integrals(eps_nodes,eps_weights,grid,order,RHO,k[1:ns,1:ns]*k[1:ns,1:ns]')

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    N = prod(length.(grid))

    for j = 1:N

        sub = ind2sub(j,Tuple(length.(grid)))

        state = Array{T,1}(undef,nx)
        for i = 1:nx
            state[i] = grid[i][sub[i]]
        end

        dr = ss_eqm.g(state)
        te = ss_eqm.h(state,zeros(ns))

        for i = 1:ny
            variables[i][j] = min(max(dr[i],lb[i]),ub[i])
        end
        for i = 1:nx
            variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
        end

    end

    if typeof(order) <: Integer
        ord = Tuple(fill(order,nx))
    else
        ord = Tuple(order)
    end

    weights = [zeros(ord.+1) for _ in 1:length(jumps_approximated)]
    scaled_weights = [zeros(ord.+1) for _ in 1:length(jumps_approximated)]

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    if node_generator == chebyshev_nodes
        cheb_weights = chebyshev_weights
    elseif node_generator == chebyshev_extrema
        cheb_weights = chebyshev_weights_extrema
    elseif node_generator == chebyshev_extended
        cheb_weights = chebyshev_weights_extended
    end

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= cheb_weights(variables[jumps_approximated[i]],grid,order,domain)
        end

        scale_chebyshev_weights!(weights,scaled_weights,integrals,jumps_approximated,ns)

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N

                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end

                init = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_chebyshev(state,scaled_weights,order,domain)
                if typeof(scheme) <: ChebyshevSchemeOBCStoch
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = ChebyshevSolutionStoch([variables[ny+1:end]; variables[1:ny]],weights,integrals,grid,order,domain,k[1:ns,:],iters,scheme.node_generator)


    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{SmolyakSchemeDet,SmolyakSchemeOBCDet})

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    layer = scheme.layer
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: SmolyakSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    T = typeof(scheme.ftol)

    grid,multi_ind = smolyak_grid(node_generator,nx,layer,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],N)
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],N)
    end

    weights  = [zeros(N) for _ in 1:length(jumps_approximated)]
    smol_iim = smolyak_inverse_interpolation_matrix(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    init = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= smolyak_weights(variables[jumps_approximated[i]],smol_iim)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_smolyak(grid[i,:],weights,multi_ind,domain)
            if typeof(scheme) <: SmolyakSchemeOBCDet
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionDet([variables[ny+1:end]; variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{SmolyakSchemeDet,SmolyakSchemeOBCDet},threads::S) where {S<:Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    layer = scheme.layer
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: SmolyakSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    T = typeof(scheme.ftol)

    grid,multi_ind = smolyak_grid(node_generator,nx,layer,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],N)
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],N)
    end

    weights  = [zeros(N) for _ in 1:length(jumps_approximated)]
    smol_iim = smolyak_inverse_interpolation_matrix_threaded(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= smolyak_weights(variables[jumps_approximated[i]],smol_iim) # threaded through the interpolation matrix
        end

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N

                init = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_smolyak(grid[i,:],weights,multi_ind,domain)
                if typeof(scheme) <: SmolyakSchemeOBCDet
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionDet([variables[ny+1:end]; variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{SmolyakSchemeStoch,SmolyakSchemeOBCStoch})

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer = scheme.layer
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: SmolyakSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    T = typeof(scheme.ftol)

    d = compute_linearization(model,initial_guess)
    GAMMA = d[1:ns,nv+1:nv+ns]
    k   = -GAMMA\d[1:ns,2*nv+1:end]
    RHO = -GAMMA\d[1:ns,1:ns]
    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),sqrt(eps()))) != ns - 1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via Smolyak methods.")
        end
    end

    grid,multi_ind = smolyak_grid(node_generator,nx,layer,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = smolyak_weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],N)
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],N)
    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    scaled_weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    smol_iim = smolyak_inverse_interpolation_matrix(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    init = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= smolyak_weights(variables[jumps_approximated[i]],smol_iim)
            scaled_weights[i] .= scale_sparse_weights(weights[i],weight_scale_factor)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_smolyak(grid[i,:],scaled_weights,multi_ind,domain)
            if typeof(scheme) <: SmolyakSchemeOBCStoch
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionStoch([variables[ny+1:end]; variables[1:ny]],weights,weight_scale_factor,grid,multi_ind,layer,domain,k,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{SmolyakSchemeStoch,SmolyakSchemeOBCStoch},threads::S) where {S<:Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer = scheme.layer
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: SmolyakSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    T = typeof(scheme.ftol)

    d = compute_linearization(model,initial_guess)
    GAMMA = d[1:ns,nv+1:nv+ns]
    k   = -GAMMA\d[1:ns,2*nv+1:end]
    RHO = -GAMMA\d[1:ns,1:ns]
    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),sqrt(eps()))) != ns - 1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via Smolyak methods.")
        end
    end

    grid,multi_ind = smolyak_grid(node_generator,nx,layer,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = smolyak_weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],N)
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],N)
    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    scaled_weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    smol_iim = smolyak_inverse_interpolation_matrix_threaded(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= smolyak_weights(variables[jumps_approximated[i]],smol_iim) # threaded through the interpolation matrix 
            scaled_weights[i] .= scale_sparse_weights(weights[i],weight_scale_factor)
        end

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N

                init = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_smolyak(grid[i,:],scaled_weights,multi_ind,domain)
                if typeof(scheme) <: SmolyakSchemeOBCStoch
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionStoch([variables[ny+1:end]; variables[1:ny]],weights,weight_scale_factor,grid,multi_ind,layer,domain,k,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{SmolyakSchemeDet,SmolyakSchemeOBCDet}) where {R<:Union{PerturbationSolutionDet,ProjectionSolutionDet}}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    layer = scheme.layer
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: SmolyakSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    ss_eqm = state_space_eqm(soln)

    T = typeof(scheme.ftol)

    grid,multi_ind = smolyak_grid(node_generator,nx,layer,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    for j = 1:N

        dr = ss_eqm.g(grid[j,:])
        te = ss_eqm.h(grid[j,:])

        for i = 1:ny
            variables[i][j] = min(max(dr[i],lb[i]),ub[i])
        end
        for i = 1:nx
            variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
        end

    end

    weights  = [zeros(N) for _ in 1:length(jumps_approximated)]
    smol_iim = smolyak_inverse_interpolation_matrix(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    init = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= smolyak_weights(variables[jumps_approximated[i]],smol_iim)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_smolyak(grid[i,:],weights,multi_ind,domain)
            if typeof(scheme) <: SmolyakSchemeOBCDet
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionDet([variables[ny+1:end]; variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{SmolyakSchemeDet,SmolyakSchemeOBCDet},threads::S) where {R<:Union{PerturbationSolutionDet,ProjectionSolutionDet},S<:Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    layer = scheme.layer
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: SmolyakSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    ss_eqm = state_space_eqm(soln)

    T = typeof(scheme.ftol)

    grid,multi_ind = smolyak_grid(node_generator,nx,layer,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    @sync Threads.@threads for t = 1:threads
        for j = t:threads:N

            dr = ss_eqm.g(grid[j,:])
            te = ss_eqm.h(grid[j,:])

            for i = 1:ny
                variables[i][j] = min(max(dr[i],lb[i]),ub[i])
            end
            for i = 1:nx
                variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
            end
    
        end

    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    smol_iim = smolyak_inverse_interpolation_matrix_threaded(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= smolyak_weights(variables[jumps_approximated[i]],smol_iim) # threaded through the interpolation matrix
        end

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N

                init = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_smolyak(grid[i,:],weights,multi_ind,domain)
                if typeof(scheme) <: SmolyakSchemeOBCDet
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionDet([variables[ny+1:end]; variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{SmolyakSchemeStoch,SmolyakSchemeOBCStoch}) where {R<:Union{PerturbationSolutionStoch,ProjectionSolutionStoch}}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer = scheme.layer
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: SmolyakSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    ss_eqm = state_space_eqm(soln)

    k = soln.k[1:ns,:]

    if typeof(soln) <: PerturbationSolution
        RHO = soln.hx[1:ns,1:ns]
    elseif typeof(soln) <: ProjectionSolution
        d = compute_linearization(model,initial_guess)
        GAMMA = d[1:ns,nv+1:nv+ns]
        k   = -GAMMA\d[1:ns,2*nv+1:end]
        RHO = -GAMMA\d[1:ns,1:ns]
    end

    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),sqrt(eps()))) != ns - 1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via Smolyak methods.")
        end
    end

    T = typeof(scheme.ftol)

    if domain == []
        if typeof(soln) <: ProjectionSolution
            domain = soln.domain
        elseif typeof(soln) <: PerturbationSolutionStoch
            soln_fo = FirstOrderSolutionStoch(soln.hbar,soln.hx,soln.k,soln.gbar,soln.gx,soln.sigma,soln.grc,soln.soln_type)
            state_vars,jump_vars = compute_variances(soln_fo)
            domain = Matrix([soln.hbar + 3*sqrt.(diag(state_vars)) soln.hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
        end
    end

    grid,multi_ind = smolyak_grid(node_generator,nx,layer,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = smolyak_weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    for j = 1:N

        dr = ss_eqm.g(grid[j,:])
        te = ss_eqm.h(grid[j,:],zeros(ns))

        for i = 1:ny
            variables[i][j] = min(max(dr[i],lb[i]),ub[i])
        end
        for i = 1:nx
            variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
        end

    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    scaled_weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    smol_iim = smolyak_inverse_interpolation_matrix(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    init = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= smolyak_weights(variables[jumps_approximated[i]],smol_iim)
            scaled_weights[i] .= scale_sparse_weights(weights[i],weight_scale_factor)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_smolyak(grid[i,:],scaled_weights,multi_ind,domain)
            if typeof(scheme) <: SmolyakSchemeOBCStoch
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionStoch([variables[ny+1:end]; variables[1:ny]],weights,weight_scale_factor,grid,multi_ind,layer,domain,k[1:ns,:],iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{SmolyakSchemeStoch,SmolyakSchemeOBCStoch},threads::S) where {R<:Union{PerturbationSolutionStoch,ProjectionSolutionStoch},S<:Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer = scheme.layer
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: SmolyakSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    ss_eqm = state_space_eqm(soln)

    k = soln.k[1:ns,:]

    if typeof(soln) <: PerturbationSolution
        RHO = soln.hx[1:ns,1:ns]
    elseif typeof(soln) <: ProjectionSolution
        d = compute_linearization(model,initial_guess)
        GAMMA = d[1:ns,nv+1:nv+ns]
        k   = -GAMMA\d[1:ns,2*nv+1:end]
        RHO = -GAMMA\d[1:ns,1:ns]
    end

    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),sqrt(eps()))) != ns - 1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via Smolyak methods.")
        end
    end

    T = typeof(scheme.ftol)

    if domain == []
        if typeof(soln) <: ProjectionSolution
            domain = soln.domain
        elseif typeof(soln) <: PerturbationSolutionStoch
            soln_fo = FirstOrderSolutionStoch(soln.hbar,soln.hx,soln.k,soln.gbar,soln.gx,soln.sigma,soln.grc,soln.soln_type)
            state_vars,jump_vars = compute_variances(soln_fo)
            domain = Matrix([soln.hbar + 3*sqrt.(diag(state_vars)) soln.hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
        end
    end

    grid,multi_ind = smolyak_grid(node_generator,nx,layer,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = smolyak_weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    @sync Threads.@threads for t = 1:threads
        for j = t:threads:N

            dr = ss_eqm.g(grid[j,:])
            te = ss_eqm.h(grid[j,:],zeros(ns))

            for i = 1:ny
                variables[i][j] = min(max(dr[i],lb[i]),ub[i])
            end
            for i = 1:nx
                variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
            end
    
        end
    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    scaled_weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    smol_iim = smolyak_inverse_interpolation_matrix_threaded(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= smolyak_weights(variables[jumps_approximated[i]],smol_iim) # threaded through the interpolation matrix
            scaled_weights[i] .= scale_sparse_weights(weights[i],weight_scale_factor)
        end

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N

                init = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_smolyak(grid[i,:],scaled_weights,multi_ind,domain)
                if typeof(scheme) <: SmolyakSchemeOBCStoch
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = SmolyakSolutionStoch([variables[ny+1:end]; variables[1:ny]],weights,weight_scale_factor,grid,multi_ind,layer,domain,k[1:ns,:],iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{PiecewiseLinearSchemeDet,PiecewiseLinearSchemeOBCDet})

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_number = scheme.node_number
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: PiecewiseLinearSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.ftol)

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

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    init = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i = 1:N
            sub = ind2sub(i,Tuple(length.(grid)))

            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_piecewise(variables,grid,state)
            if typeof(scheme) <: PiecewiseLinearSchemeOBCDet
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionDet([variables[ny+1:end]; variables[1:ny]],grid,domain,iters)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{PiecewiseLinearSchemeDet,PiecewiseLinearSchemeOBCDet},threads::S) where {S<:Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_number = scheme.node_number
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: PiecewiseLinearSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.ftol)

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

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N
                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_piecewise(variables,grid,state)
                if typeof(scheme) <: PiecewiseLinearSchemeOBCDet
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionDet([variables[ny+1:end]; variables[1:ny]],grid,domain,iters)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{PiecewiseLinearSchemeStoch,PiecewiseLinearSchemeOBCStoch})

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_number = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: PiecewiseLinearSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.ftol)

    d = compute_linearization(model,initial_guess)
    GAMMA = d[1:ns,nv+1:nv+ns]
    k   = -GAMMA\d[1:ns,2*nv+1:end]
    RHO = -GAMMA\d[1:ns,1:ns]
    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
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

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    init = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i = 1:N
            sub = ind2sub(i,Tuple(length.(grid)))

            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_piecewise(variables,grid,state,integrals)
            if typeof(scheme) <: PiecewiseLinearSchemeOBCStoch
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionStoch([variables[ny+1:end]; variables[1:ny]],grid,domain,k,iters)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{PiecewiseLinearSchemeStoch,PiecewiseLinearSchemeOBCStoch},threads::S) where {S<:Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_number = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: PiecewiseLinearSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    T = typeof(scheme.ftol)

    d = compute_linearization(model,initial_guess)
    GAMMA = d[1:ns,nv+1:nv+ns]
    k   = -GAMMA\d[1:ns,2*nv+1:end]
    RHO = -GAMMA\d[1:ns,1:ns]
    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
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

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N
                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_piecewise(variables,grid,state,integrals)
                if typeof(scheme) <: PiecewiseLinearSchemeOBCStoch
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end
                
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionStoch([variables[ny+1:end]; variables[1:ny]],grid,domain,k,iters)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{PiecewiseLinearSchemeDet,PiecewiseLinearSchemeOBCDet}) where {R<:Union{PerturbationSolutionDet,ProjectionSolutionDet}}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    initial_guess = scheme.initial_guess
    node_number = scheme.node_number
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: PiecewiseLinearSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    ss_eqm = state_space_eqm(soln)

    T = typeof(scheme.ftol)

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

        dr = ss_eqm.g(state)
        te = ss_eqm.h(state)

        for i = 1:ny
            variables[i][j] = min(max(dr[i],lb[i]),ub[i])
        end
        for i = 1:nx
            variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
        end

    end

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    init = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i = 1:N
            sub = ind2sub(i,Tuple(length.(grid)))

            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_piecewise(variables,grid,state)
            if typeof(scheme) <: PiecewiseLinearSchemeOBCDet
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionDet([variables[ny+1:end]; variables[1:ny]],grid,domain,iters)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{PiecewiseLinearSchemeDet,PiecewiseLinearSchemeOBCDet},threads::S) where {R<:Union{PerturbationSolutionDet,ProjectionSolutionDet},S<:Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    initial_guess = scheme.initial_guess
    node_number = scheme.node_number
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: PiecewiseLinearSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    ss_eqm = state_space_eqm(soln)

    T = typeof(scheme.ftol)

    grid = Array{Array{T,1},1}(undef,nx)
    for i = 1:nx
        grid[i] = piecewise_linear_nodes(node_number[i],domain[:,i])
    end

    N = prod(length.(grid))

    variables = Array{Array{T,nx},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    @sync Threads.@threads for t = 1:threads
        for j = t:threads:N
            sub = ind2sub(j,Tuple(length.(grid)))

            state = Array{T,1}(undef,nx)
            for i = 1:nx
                state[i] = grid[i][sub[i]]
            end

            dr = ss_eqm.g(state)
            te = ss_eqm.h(state)

            for i = 1:ny
                variables[i][j] = min(max(dr[i],lb[i]),ub[i])
            end
            for i = 1:nx
                variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
            end
    
        end
    end

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N
                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_piecewise(variables,grid,state)
                if typeof(scheme) <: PiecewiseLinearSchemeOBCDet
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end
    
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionDet([variables[ny+1:end]; variables[1:ny]],grid,domain,iters)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{PiecewiseLinearSchemeStoch,PiecewiseLinearSchemeOBCStoch}) where {R<:Union{PerturbationSolutionStoch,ProjectionSolutionStoch}}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    initial_guess = scheme.initial_guess
    node_number = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: PiecewiseLinearSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    ss_eqm = state_space_eqm(soln)

    k = soln.k[1:ns,:]

    if typeof(soln) <: PerturbationSolution
        RHO = soln.hx[1:ns,1:ns]
    elseif typeof(soln) <: ProjectionSolution
        d = compute_linearization(model,initial_guess)
        GAMMA = d[1:ns,nv+1:nv+ns]
        k   = -GAMMA\d[1:ns,2*nv+1:end]
        RHO = -GAMMA\d[1:ns,1:ns]
    end

    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),sqrt(eps()))) != ns - 1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via Smolyak methods.")
        end
    end

    T = typeof(scheme.ftol)

    if domain == []
        if typeof(soln) <: ProjectionSolution
            domain = soln.domain
        elseif typeof(soln) <: PerturbationSolutionStoch
            soln_fo = FirstOrderSolutionStoch(soln.hbar,soln.hx,soln.k,soln.gbar,soln.gx,soln.sigma,soln.grc,soln.soln_type)
            state_vars,jump_vars = compute_variances(soln_fo)
            domain = Matrix([soln.hbar + 3*sqrt.(diag(state_vars)) soln.hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
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
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    for j = 1:N
        sub = ind2sub(j,Tuple(length.(grid)))

        for i = 1:nx
            state[i] = grid[i][sub[i]]
        end

        dr = ss_eqm.g(state)
        te = ss_eqm.h(state,zeros(ns))

        for i = 1:ny
            variables[i][j] = min(max(dr[i],lb[i]),ub[i])
        end
        for i = 1:nx
            variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
        end

    end

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    init = Array{T,1}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i = 1:N
            sub = ind2sub(i,Tuple(length.(grid)))

            for j = 1:nx
                state[j] = grid[j][sub[j]]
            end
            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_piecewise(variables,grid,state,integrals)
            if typeof(scheme) <: PiecewiseLinearSchemeOBCStoch
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end
        
            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionStoch([variables[ny+1:end]; variables[1:ny]],grid,domain,k[1:ns,:],iters)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{PiecewiseLinearSchemeStoch,PiecewiseLinearSchemeOBCStoch},threads::S) where {R<:Union{PerturbationSolutionStoch,ProjectionSolutionStoch},S<:Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    initial_guess = scheme.initial_guess
    node_number = scheme.node_number
    num_quad_nodes = scheme.num_quad_nodes
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: PiecewiseLinearSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    if length(node_number) != nx
        error("The number of nodes is needed for each state and only each state variable.")
    end

    ss_eqm = state_space_eqm(soln)

    k = soln.k[1:ns,:]

    if typeof(soln) <: PerturbationSolution
        RHO = soln.hx[1:ns,1:ns]
    elseif typeof(soln) <: ProjectionSolution
        d = compute_linearization(model,initial_guess)
        GAMMA = d[1:ns,nv+1:nv+ns]
        k   = -GAMMA\d[1:ns,2*nv+1:end]
        RHO = -GAMMA\d[1:ns,1:ns]
    end

    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),sqrt(eps()))) != ns - 1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via Smolyak methods.")
        end
    end

    T = typeof(scheme.ftol)

    if domain == []
        if typeof(soln) <: ProjectionSolution
            domain = soln.domain
        elseif typeof(soln) <: PerturbationSolutionStoch
            soln_fo = FirstOrderSolutionStoch(soln.hbar,soln.hx,soln.k,soln.gbar,soln.gx,soln.sigma,soln.grc,soln.soln_type)
            state_vars,jump_vars = compute_variances(soln_fo)
            domain = Matrix([soln.hbar + 3*sqrt.(diag(state_vars)) soln.hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
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
    for i = 1:nv
        variables[i] = zeros(Tuple(length.(grid)))
    end

    @sync Threads.@threads for t = 1:threads
        for j = t:threads:N
            sub = ind2sub(j,Tuple(length.(grid)))

            state = Array{T,1}(undef,nx)
            for i = 1:nx
                state[i] = grid[i][sub[i]]
            end

            dr = ss_eqm.g(state)
            te = ss_eqm.h(state,zeros(ns))

            for i = 1:ny
                variables[i][j] = min(max(dr[i],lb[i]),ub[i])
            end
            for i = 1:nx
                variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
            end
    
        end
    end

    new_variables = [zeros(Tuple(length.(grid))) for _ in 1:nv]

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N
                sub = ind2sub(i,Tuple(length.(grid)))

                state = Array{T,1}(undef,nx)
                for j = 1:nx
                    state[j] = grid[j][sub[j]]
                end
                init = Array{T,1}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_piecewise(variables,grid,state,integrals)
                if typeof(scheme) <: PiecewiseLinearSchemeOBCStoch
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end
    
                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = PiecewiseLinearSolutionStoch([variables[ny+1:end]; variables[1:ny]],grid,domain,k[1:ns,:],iters)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{HyperbolicCrossSchemeDet,HyperbolicCrossSchemeOBCDet})

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    layer = scheme.layer
    n = scheme.n
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: HyperbolicCrossSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    T = typeof(scheme.ftol)

    grid,multi_ind = hyperbolic_cross_grid(node_generator,nx,layer,n,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],N)
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],N)
    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    init = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= hyperbolic_cross_weights(variables[jumps_approximated[i]],hcross_iim)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_hcross(grid[i,:],weights,multi_ind,domain)
            if typeof(scheme) <: HyperbolicCrossSchemeOBCDet
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = HyperbolicCrossSolutionDet([variables[ny+1:end]; variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{HyperbolicCrossSchemeDet,HyperbolicCrossSchemeOBCDet},threads::S) where {S<:Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    layer = scheme.layer
    n = scheme.n
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: HyperbolicCrossSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    T = typeof(scheme.ftol)

    grid,multi_ind = hyperbolic_cross_grid(node_generator,nx,layer,n,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],N)
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],N)
    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix_threaded(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= hyperbolic_cross_weights(variables[jumps_approximated[i]],hcross_iim) # threaded through the interpolation matrix
        end

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N

                init = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_hcross(grid[i,:],weights,multi_ind,domain)
                if typeof(scheme) <: HyperbolicCrossSchemeOBCDet
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = HyperbolicCrossSolutionDet([variables[ny+1:end]; variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{HyperbolicCrossSchemeStoch,HyperbolicCrossSchemeOBCStoch})

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer = scheme.layer
    n = scheme.n
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: HyperbolicCrossSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    T = typeof(scheme.ftol)

    d = compute_linearization(model,initial_guess)
    GAMMA = d[1:ns,nv+1:nv+ns]
    k   = -GAMMA\d[1:ns,2*nv+1:end]
    RHO = -GAMMA\d[1:ns,1:ns]
    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),sqrt(eps()))) != ns - 1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via Hyperbolic cross methods.")
        end
    end

    grid,multi_ind = hyperbolic_cross_grid(node_generator,nx,layer,n,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = hcross_weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],N)
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],N)
    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    scaled_weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    init = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= hyperbolic_cross_weights(variables[jumps_approximated[i]],hcross_iim)
            scaled_weights[i] .= scale_sparse_weights(weights[i],weight_scale_factor)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_hcross(grid[i,:],scaled_weights,multi_ind,domain)
            if typeof(scheme) <: HyperbolicCrossSchemeOBCStoch
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = HyperbolicCrossSolutionStoch([variables[ny+1:end]; variables[1:ny]],weights,weight_scale_factor,grid,multi_ind,layer,domain,k,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,scheme::Union{HyperbolicCrossSchemeStoch,HyperbolicCrossSchemeOBCStoch},threads::S) where {S<:Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer = scheme.layer
    n = scheme.n
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: HyperbolicCrossSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    T = typeof(scheme.ftol)

    d = compute_linearization(model,initial_guess)
    GAMMA = d[1:ns,nv+1:nv+ns]
    k   = -GAMMA\d[1:ns,2*nv+1:end]
    RHO = -GAMMA\d[1:ns,1:ns]
    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),sqrt(eps()))) != ns - 1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via Hyperbolic cross methods.")
        end
    end

    grid,multi_ind = hyperbolic_cross_grid(node_generator,nx,layer,n,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = hcross_weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:ny
        variables[i] = fill(initial_guess[nx+i],N)
    end
    for i = 1:nx
        variables[ny+i] = fill(initial_guess[i],N)
    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    scaled_weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix_threaded(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= hyperbolic_cross_weights(variables[jumps_approximated[i]],hcross_iim) # threaded through the interpolation matrix
            scaled_weights[i] .= scale_sparse_weights(weights[i],weight_scale_factor)
        end

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N

                init = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_hcross(grid[i,:],scaled_weights,multi_ind,domain)
                if typeof(scheme) <: HyperbolicCrossSchemeOBCStoch
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = HyperbolicCrossSolutionStoch([variables[ny+1:end]; variables[1:ny]],weights,weight_scale_factor,grid,multi_ind,layer,domain,k,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{HyperbolicCrossSchemeDet,HyperbolicCrossSchemeOBCDet}) where {R<:Union{PerturbationSolutionDet,ProjectionSolutionDet}}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    layer = scheme.layer
    n = scheme.n
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: HyperbolicCrossSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    ss_eqm = state_space_eqm(soln)

    T = typeof(scheme.ftol)

    grid,multi_ind = hyperbolic_cross_grid(node_generator,nx,layer,n,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    for j = 1:N

        dr = ss_eqm.g(grid[j,:])
        te = ss_eqm.h(grid[j,:])

        for i = 1:ny
            variables[i][j] = min(max(dr[i],lb[i]),ub[i])
        end
        for i = 1:nx
            variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
        end

    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    init = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= hyperbolic_cross_weights(variables[jumps_approximated[i]],hcross_iim)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_hcross(grid[i,:],weights,multi_ind,domain)
            if typeof(scheme) <: HyperbolicCrossSchemeOBCDet
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = HyperbolicCrossSolutionDet([variables[ny+1:end]; variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{HyperbolicCrossSchemeDet,HyperbolicCrossSchemeOBCDet},threads::S) where {R<:Union{PerturbationSolutionDet,ProjectionSolutionDet},S<:Integer}

    nx = model.number_states
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    layer = scheme.layer
    n = scheme.n
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: HyperbolicCrossSchemeOBCDet
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    ss_eqm = state_space_eqm(soln)

    T = typeof(scheme.ftol)

    grid,multi_ind = hyperbolic_cross_grid(node_generator,nx,layer,n,domain)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    @sync Threads.@threads for t = 1:threads
        for j = t:threads:N

            dr = ss_eqm.g(grid[j,:])
            te = ss_eqm.h(grid[j,:])

            for i = 1:ny
                variables[i][j] = min(max(dr[i],lb[i]),ub[i])
            end
            for i = 1:nx
                variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
            end
    
        end

    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix_threaded(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= hyperbolic_cross_weights(variables[jumps_approximated[i]],hcross_iim) # threaded through the interpolation matrix
        end

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N

                init = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_hcross(grid[i,:],weights,multi_ind,domain)
                if typeof(scheme) <: HyperbolicCrossSchemeOBCDet
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = HyperbolicCrossSolutionDet([variables[ny+1:end]; variables[1:ny]],weights,grid,multi_ind,layer,domain,iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{HyperbolicCrossSchemeStoch,HyperbolicCrossSchemeOBCStoch}) where {R<:Union{PerturbationSolutionStoch,ProjectionSolutionStoch}}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer = scheme.layer
    n = scheme.n
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: HyperbolicCrossSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    ss_eqm = state_space_eqm(soln)

    k = soln.k[1:ns,:]

    if typeof(soln) <: PerturbationSolution
        RHO = soln.hx[1:ns,1:ns]
    elseif typeof(soln) <: ProjectionSolution
        d = compute_linearization(model,initial_guess)
        GAMMA = d[1:ns,nv+1:nv+ns]
        k   = -GAMMA\d[1:ns,2*nv+1:end]
        RHO = -GAMMA\d[1:ns,1:ns]
    end

    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),sqrt(eps()))) != ns - 1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via Hyperbolic cross methods.")
        end
    end

    T = typeof(scheme.ftol)

    if domain == []
        if typeof(soln) <: ProjectionSolution
            domain = soln.domain
        elseif typeof(soln) <: PerturbationSolutionStoch
            soln_fo = FirstOrderSolutionStoch(soln.hbar,soln.hx,soln.k,soln.gbar,soln.gx,soln.sigma,soln.grc,soln.soln_type)
            state_vars,jump_vars = compute_variances(soln_fo)
            domain = Matrix([soln.hbar + 3*sqrt.(diag(state_vars)) soln.hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
        end
    end

    grid,multi_ind = hyperbolic_cross_grid(node_generator,nx,layer,n,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = hcross_weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    for j = 1:N

        dr = ss_eqm.g(grid[j,:])
        te = ss_eqm.h(grid[j,:],zeros(ns))

        for i = 1:ny
            variables[i][j] = min(max(dr[i],lb[i]),ub[i])
        end
        for i = 1:nx
            variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
        end

    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    scaled_weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    init = Array{T}(undef,nv)

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= hyperbolic_cross_weights(variables[jumps_approximated[i]],hcross_iim)
            scaled_weights[i] .= scale_sparse_weights(weights[i],weight_scale_factor)
        end

        for i = 1:N

            for j = 1:nv
                init[j] = variables[j][i]
            end

            projection_equations = model.closure_function_hcross(grid[i,:],scaled_weights,multi_ind,domain)
            if typeof(scheme) <: HyperbolicCrossSchemeOBCStoch
                nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
            else
                nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
            end

            for j = 1:nv
                new_variables[j][i] = nlsoln.zero[j]
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = HyperbolicCrossSolutionStoch([variables[ny+1:end]; variables[1:ny]],weights,weight_scale_factor,grid,multi_ind,layer,domain,k[1:ns,:],iters,scheme.node_generator)

    return soln

end

function solve_nonlinear(model::REModel,soln::R,scheme::Union{HyperbolicCrossSchemeStoch,HyperbolicCrossSchemeOBCStoch},threads::S) where {R<:Union{PerturbationSolutionStoch,ProjectionSolutionStoch},S<:Integer}

    nx = model.number_states
    ns = model.number_shocks
    ny = model.number_jumps
    nv = nx + ny

    jumps_approximated = model.jumps_approximated

    initial_guess = scheme.initial_guess
    node_generator = scheme.node_generator
    num_quad_nodes = scheme.num_quad_nodes
    layer = scheme.layer
    n = scheme.n
    domain = scheme.domain

    ub = [Inf for _ in 1:nv]
    lb = [-Inf for _ in 1:nv]
    if typeof(scheme) <: HyperbolicCrossSchemeOBCStoch
        ub[1:ny]    .= scheme.ub[nx+1:nv]
        ub[ny+1:nv] .= scheme.ub[1:nx]
        lb[1:ny]    .= scheme.lb[nx+1:nv]
        lb[ny+1:nv] .= scheme.lb[1:nx]
    end

    ss_eqm = state_space_eqm(soln)

    k = soln.k[1:ns,:]

    if typeof(soln) <: PerturbationSolution
        RHO = soln.hx[1:ns,1:ns]
    elseif typeof(soln) <: ProjectionSolution
        d = compute_linearization(model,initial_guess)
        GAMMA = d[1:ns,nv+1:nv+ns]
        k   = -GAMMA\d[1:ns,2*nv+1:end]
        RHO = -GAMMA\d[1:ns,1:ns]
    end

    if !isdiag(RHO .> sqrt(eps()))
        error("This solver requires the shocks to be AR(1) processes.")
    end
    for i = 1:ns
        if sum(isless.(abs.(k[i,:]),sqrt(eps()))) != ns - 1 # Make sure there is only one shock per equation
            error("Models with correlated shocks cannot yet be solved via Hyperbolic cross methods.")
        end
    end

    T = typeof(scheme.ftol)

    if domain == []
        if typeof(soln) <: ProjectionSolution
            domain = soln.domain
        elseif typeof(soln) <: PerturbationSolutionStoch
            soln_fo = FirstOrderSolutionStoch(soln.hbar,soln.hx,soln.k,soln.gbar,soln.gx,soln.sigma,soln.grc,soln.soln_type)
            state_vars,jump_vars = compute_variances(soln_fo)
            domain = Matrix([soln.hbar + 3*sqrt.(diag(state_vars)) soln.hbar - 3*sqrt.(diag(state_vars))]')   # dimension are 2*nx
        end
    end

    grid,multi_ind = hyperbolic_cross_grid(node_generator,nx,layer,n,domain)
    (eps_nodes,eps_weights) = hermite(num_quad_nodes)
    weight_scale_factor = hcross_weight_scale_factors(eps_nodes,eps_weights,multi_ind,nx,grid,RHO,k)

    N = size(grid,1)

    variables = Array{Array{T,1},1}(undef,nv)
    for i = 1:nv
        variables[i] = zeros(N)
    end

    @sync Threads.@threads for t = 1:threads
        for j = t:threads:N

            dr = ss_eqm.g(grid[j,:])
            te = ss_eqm.h(grid[j,:],zeros(ns))

            for i = 1:ny
                variables[i][j] = min(max(dr[i],lb[i]),ub[i])
            end
            for i = 1:nx
                variables[ny+i][j] = min(max(te[i],lb[ny+i]),ub[ny+i])
            end
    
        end
    end

    weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    scaled_weights = [zeros(N) for _ in 1:length(jumps_approximated)]
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix_threaded(grid,multi_ind,domain)

    new_variables = [zeros(N) for _ in 1:nv]

    iters = 0
    len = Inf
    while len > scheme.xtol && iters <= scheme.maxiters

        for i in eachindex(jumps_approximated)
            weights[i] .= hyperbolic_cross_weights(variables[jumps_approximated[i]],hcross_iim) # threaded through the interpolation matrix
            scaled_weights[i] .= scale_sparse_weights(weights[i],weight_scale_factor)
        end

        @sync Threads.@threads for t = 1:threads
            for i = t:threads:N

                init = Array{T}(undef,nv)
                for j = 1:nv
                    init[j] = variables[j][i]
                end

                projection_equations = model.closure_function_hcross(grid[i,:],scaled_weights,multi_ind,domain)
                if typeof(scheme) <: HyperbolicCrossSchemeOBCStoch
                    nlsoln = nlboxsolve(projection_equations,init,lb,ub,xtol = scheme.xtol,ftol = scheme.ftol,iterations = scheme.maxiters,method = scheme.method)
                else
                    nlsoln = nlsolve(projection_equations,init,xtol = scheme.ftol,iterations = scheme.maxiters,inplace = :true,method = scheme.method)
                end

                for j = 1:nv
                    new_variables[j][i] = nlsoln.zero[j]
                end
            end
        end

        len = zero(T)
        for j = 1:nv
            len = mylen(len,new_variables[j],variables[j])
        end

        for j = 1:nv
            variables[j] .= new_variables[j]
        end

        iters += 1

    end

    soln = HyperbolicCrossSolutionStoch([variables[ny+1:end]; variables[1:ny]],weights,weight_scale_factor,grid,multi_ind,layer,domain,k[1:ns,:],iters,scheme.node_generator)

    return soln

end

"""
Solves a model given a solution scheme.

Exported function.
"""
function solve_model(model::REModel,scheme::PerturbationScheme)

    if scheme.order == "first" && (model.solvers in ("Any", "Linear", "Perturbation"))
        soln = solve_first_order(model, scheme)
        return soln
    elseif scheme.order == "second" && (model.solvers in ("Any", "Perturbation"))
        soln = solve_second_order(model, scheme)
        return soln
    elseif scheme.order == "third" && (model.solvers in ("Any", "Perturbation"))
        soln = solve_third_order(model, scheme)
        return soln
    elseif scheme.order == "fourth" && (model.solvers in ("Any", "Perturbation"))
        soln = solve_fourth_order(model, scheme)
        return soln
    else
        error("The chosen order has not been implemented or the solution scheme conflicts with the solvers specified in the model file.")
    end

end

function solve_model(model::REModel,scheme::P) where {P<:ProjectionScheme}

    if model.solvers in ("Any", "Projection")
        if model.number_shocks != 0 && typeof(scheme) <: Union{ProjectionSchemeDet,ProjectionSchemeOBCDet}
            error("Stochastic model, but deterministic SolutionScheme.")
        elseif model.number_shocks == 0 && typeof(scheme) <: Union{ProjectionSchemeStoch,ProjectionSchemeOBCStoch}
            error("Deterministic model, but stochastic SolutionScheme.")
        end

        soln = solve_nonlinear(model,scheme)
        return soln
    else
        error("The solution scheme conflicts with the solvers specified in the model file.")
    end

end

function solve_model(model::REModel,scheme::P,threads::S) where {P<:ProjectionScheme,S<:Integer}

    if model.solvers in ("Any", "Projection")
        if model.number_shocks != 0 && typeof(scheme) <: Union{ProjectionSchemeDet,ProjectionSchemeOBCDet}
            error("Stochastic model, but deterministic SolutionScheme.")
        elseif model.number_shocks == 0 && typeof(scheme) <: Union{ProjectionSchemeStoch,ProjectionSchemeOBCStoch}
            error("Deterministic model, but stochastic SolutionScheme.")
        end

        if threads < 0
            error("Number of threads cannot be negative.")
        elseif threads == 0
            soln = solve_nonlinear(model,scheme)
            return soln
        else
            soln = solve_nonlinear(model,scheme,threads)
            return soln
        end
    else
        error("The solution scheme conflicts with the solvers specified in the model file.")
    end

end

function solve_model(model::REModel,soln::ModelSolution,scheme::P) where {P<:ProjectionScheme}

    if model.solvers in ("Any", "Projection")
        if model.number_shocks != 0 && typeof(scheme) <: Union{ProjectionSchemeDet,ProjectionSchemeOBCDet}
            error("Stochastic model, but deterministic SolutionScheme.")
        elseif model.number_shocks == 0 && typeof(scheme) <: Union{ProjectionSchemeStoch,ProjectionSchemeOBCStoch}
            error("Deterministic model, but stochastic SolutionScheme.")
        end

        soln = solve_nonlinear(model,soln,scheme)
        return soln
    else
        error("The solution scheme conflicts with the solvers specified in the model file.")
    end

end

function solve_model(model::REModel,soln::ModelSolution,scheme::P,threads::S) where {P<:ProjectionScheme,S<:Integer}

    if model.solvers in ("Any", "Projection")
        if model.number_shocks != 0 && typeof(scheme) <: Union{ProjectionSchemeDet,ProjectionSchemeOBCDet}
            error("Stochastic model, but deterministic SolutionScheme.")
        elseif model.number_shocks == 0 && typeof(scheme) <: Union{ProjectionSchemeStoch,ProjectionSchemeOBCStoch}
            error("Deterministic model, but stochastic SolutionScheme.")
        end

        if threads < 0
            error("Number of threads cannot be negative.")
        elseif threads == 0
            soln = solve_nonlinear(model,soln,scheme)
            return soln
        else
            soln = solve_nonlinear(model,soln,scheme,threads)
            return soln
        end
    else
        error("The solution scheme conflicts with the solvers specified in the model file.")
    end

end