function eye(tpe::Type,n::S) where {S <: Integer}

    if tpe <: Number
        return Matrix{tpe}(I,n,n)
    end

end

eye(n::Integer) = eye(Float64,n::Integer)

function tracem(x::Array{T,2}) where {T <: Real}

  # Computes the matrix-trace as defined by Gomme and Klein (2011)

  # We require the number of rows to be greater than the number of columns, so
  # that m is greater than one.

  trans = false

  if size(x,1) < size(x,2)
    x = Matrix(x')
    trans = true
  end

  n = size(x,2)
  m = Int(size(x,1)/n)

  y = zeros(m,1)  # We want this to be a 2-d array for subsequent matrix multiplication

  @inbounds for i = 1:m
    @views y[i,1] = tr(x[(n*(i-1)+1):i*n,1:n])
  end

  if trans == true
    y = y'
  end

  return Matrix(y)

end

function kron_prod_times_vector(a::AbstractArray{T,2},b::AbstractArray{T,2},v::AbstractArray{T,1}) where {T <: Real}

    # This function efficiently computes kron(a,b)*v where a and b are matrices
    # and v is a vector.

    (n1,n2) = size(a)
    (n3,n4) = size(b)
    n5      = length(v)

    if n5 != n2*n4
        error("a, b, and v do not have the correct sizes")
    end

    p = vec(b*reshape(v,n4,n2)*a')

    return p

end

function vector_times_kron_prod(v::AbstractArray{T,1},a::AbstractArray{T,2},b::AbstractArray{T,2}) where {T <: Real}

    product = Matrix(kron_prod_times_vector(a',b',v')')

    return product

end

function kron_prod_times_matrix(a::AbstractArray{T,2},b::AbstractArray{T,2},v::AbstractArray{T,2}) where {T <: Real}

    # This function efficiently computes kron(a,b)*v where a, b, and v are
    # conformable matrices.

    (n1,n2) = size(a)
    (n3,n4) = size(b)
    (n5,n6) = size(v)

    if n5 != n2*n4
        error("a, b, and v do not have the correct sizes")
    end

    p = zeros(n1*n3,n6)
    for i = 1:n6
        @views p[:,i] = kron_prod_times_vector(a,b,v[:,i])
    end

    return p

end

function matrix_times_kron_prod(v::AbstractArray{T,2},a::AbstractArray{T,2},b::AbstractArray{T,2}) where {T <: Real}

    product = Matrix(kron_prod_times_matrix(a',b',v')')

    return product

end

function kron_prod_times_vector(a::AbstractArray{T,2},b::AbstractArray{T,2},c::AbstractArray{T,2},v::AbstractArray{T,1}) where {T <: Real}

    # This function efficiently computes kron(a,kron(b,c))*v where a, b, and c
    # are matrices and v is a vector.  Included for completeness, but not
    # actually used.

    (n1,n2) = size(a)
    (n3,n4) = size(b)
    (n5,n6) = size(c)
    n7      = length(v)

    if n7 != n2*n4*n6
        error("a, b, c, and v do not have the correct sizes")
    end

    v_tilda = reshape(v,n4*n6,n2)*a'
    p = vec(kron_prod_times_matrix(b,c,v_tilda))

    return p

end

function vector_times_kron_prod(v::AbstractArray{T,1},a::AbstractArray{T,2},b::AbstractArray{T,2},c::AbstractArray{T,2}) where {T <: Real}

    product = Matrix(kron_prod_times_vector(a',b',c',v')')

    return product

end

function kron_prod_times_matrix(a::AbstractArray{T,2},b::AbstractArray{T,2},c::AbstractArray{T,2},v::AbstractArray{T,2}) where {T <: Real}

    # This function efficiently computes kron(a,kron(b,c))*v where a, b, c and
    # v are comformable matrices.

    (n1,n2) = size(a)
    (n3,n4) = size(b)
    (n5,n6) = size(c)
    (n7,n8) = size(v)

    if n7 != n2*n4*n6
        error("a, b, c, and v do not have the correct sizes")
    end

    p = zeros(n1*n3*n5,n8)
    for i = 1:n8
        @views v_tilda = reshape(v[:,i],n4*n6,n2)*a'
        p[:,i] .= vec(kron_prod_times_matrix(b,c,v_tilda))
    end

    return p

end

function matrix_times_kron_prod(v::AbstractArray{T,2},a::AbstractArray{T,2},b::AbstractArray{T,2},c::AbstractArray{T,2}) where {T <: Real}

    product = Matrix(kron_prod_times_matrix(a',b',c',v')')

    return product

end

function dsylvester(a::AbstractArray{T,2}, b::AbstractArray{T,2}, c::Union{AbstractArray{T,1},AbstractArray{T,2}}) where {T <: Real}

    #= Uses the Hessenberg-Schur method to find the bounded solution of the
       discrete Sylvester equation:

             X + A*X*B = C

       Based on Golub, Nash, and Van Loan (1979).

    =#

    n = size(a, 1)
    m = size(b, 1)
    x = zeros(size(c))

    (s,u) = schur(Matrix(b'))
    (v,t) = hessenberg(a)

    c = v'*c*u

    j = m
    while j > 0
        j1 = j
        if j == 1
            block_size = 1
        elseif isequal(s[j,j-1],0.0) == false
            block_size = 2
            j -= 1
        else;
            block_size = 1
        end
        @views ajj = kron(s[j:j1,j:j1],t) + I
        @views rhs = vec(c[:,j:j1])
        if j1 < m
            @views rhs2 = t*(x[:,(j+1):m]*s[j:j1,(j+1):m]')
            rhs -= vec(rhs2)
        end
        w = ajj\rhs
        @views x[:,j] = w[1:n]
        if block_size == 2
            @views x[:,j1] = w[(n+1):2*n]
        end
        j -= 1
    end

    x = v*x*u'

    return x

end

function trm(x::AbstractArray{T,2}) where {T <: Real}

    # Computes the matrix trace as defined by Binning (2013).  Used for
    # computing the second-order terms, hss, gss.

    (n1,n2) = size(x)
    k = Int(round(sqrt(n2)))

    y = zeros(n1)
    for i = 1:k
        @views y += x[:,i+(i-1)*k]
    end

    return y

end

function trm2(x::AbstractArray{T,2}) where {T <: Real}

    # Computes the matrix trace as defined by Binning (2013).  Used for
    # computing the third-order terms, hssx, gssx.

    (n1,n2) = size(x)
    k = Int(round(n2^(1//3)))

    y = zeros(n1,k)
    for j = 1:k
        for i = 1:k
            @views y[:,j] += x[:,(j-1)+i+(i-1)*k^2]
        end
    end

    return y

end

function create_omega3(n::S) where {S <: Integer}

    # Creates the combination matrix for a third-order perturbation as defined in Levintal (2017).

    # This function is a simplified version of the create_OMEGA function originally written in
    # Matlab by Oren Levintal for his paper "Fifth Order Perturbation Solution to DSGE Models"
    # published in the Journal of Economic Dynamics and Control, 2017.  Permission to translate
    # this function into Julia and release it within the SolveDSGE module was granted by Oren
    # Levintal on February 5, 2020.

    ind = [1:n^3;]
    M = reshape(ind,1,n,n,n)
    Ix = eye(S,n^3)
    omega3 = (reshape(Ix[:,PermutedDimsArray(M,[1,4,2,3])],n^3,n^3)
           + reshape(Ix[:,PermutedDimsArray(M,[1,2,4,3])],n^3,n^3)
           + reshape(Ix[:,PermutedDimsArray(M,[1,2,3,4])],n^3,n^3))

    return omega3

end

function kron_prod_times_vector(A::Union{Array{Array{T,2},1},Array{Array{Complex{T},2},1}},x::Union{Array{T,1},Array{Complex{T},1}},n::Array{S,1},p::S) where {T <: Real, S <: Integer}

    # Computes y = (A[p] * A[p-1] * ... * A[1] )*x

    N = prod(n[1:p])
    z = copy(x)
    for i = 1:p
        z = (A[i]*reshape(z,n[i],Int(N/n[i])))'
    end
    y = reshape(Matrix(z),N)

    return y

end

function kron_prod_times_matrix(A::Union{Array{Array{T,2},1},Array{Array{Complex{T},2},1}},x::Union{Array{T,2},Array{Complex{T},2}},n::Array{S,1},p::S) where {T <: Real, S <: Integer}

    # Computes y = (A[p] * A[p-1] * ... * A[1] )*x

    N = prod(n[1:p])
    y = Array{Complex{T}}(undef,N,size(x,2))
    for j = 1:size(x,2)
        @views z = x[:,j]
        for i = 1:p
            z = (A[i]*reshape(z,n[i],Int(N/n[i])))'
        end
        y[:,j] .= reshape(Matrix(z),N)
    end

    return y

end

function KPShiftSolve(TT::Union{Array{Array{T,2},1},Array{Array{Complex{T},2},1}},n::Array{S,1},c::Union{Array{T,1},Array{Complex{T},1}},lambda::T,alpha::Union{T,Complex{T}}) where {T <: Real, S <: Integer}

    p = length(n)
    N = prod(n)

    c = copy(c)

    TT[p] = alpha*TT[p]
    if p == 1
        y = (TT[1] + lambda*Matrix{T}(I,n[1],n[1]))\c
    else
        y = Array{Complex{T}}(undef,N)
        mp = Int(N/n[p])
        for i = n[p]:-1:1
            idx = ((i-1)*mp+1):(i*mp)
            y[idx] = KPShiftSolve(TT[1:(p-1)],n[1:(p-1)],c[idx],lambda,TT[p][i,i])
            z = kron_prod_times_vector(TT[1:p],y[idx],n,p-1)
            for j = 1:(i-1)
                jdx = ((j-1)*mp+1):(j*mp)
                c[jdx] = c[jdx] - TT[p][j,i]*z
            end
        end
    end

    return y

end

function martin_van_loan(a::Array{T,2},b::Array{T,2},c::Array{T,2},d::Array{T,2},k::S) where {T <: Real, S <: Integer}

    #= Uses a recursive Schur method to find the bounded solution of the Sylvester
       equation:

       AX + B*X*(Kron^(k)C) = D

       Based on Martin and Van Loan (2006).  This is a simplified implementation
       of their algorithm, but it captures most of the gains over Golub, Nash,
       and van Loan (1979).

    =#

    a = copy(a)
    b = a\copy(b)
    c = copy(c)
    d = a\copy(d)

    (v,s) = hessenberg(b)       # v*s*v' = b
    (t,q) = schur(complex(c'))  # q*t*q' = c'

    v = Matrix(v)

    p = k + 2
    TT = Array{typeof(t)}(undef,p)
    TT[1] = s
    for i = 2:p
        TT[i] = conj(t)
    end

    Q     = fill(q,k+1)
    inv_Q = fill(Matrix(q'),k+1)

    n = fill(size(c,1),p)
    n[1] = size(b,1)
    N = prod(n)

    lambda = 1.0

    e = vec(v'*kron_prod_times_matrix(inv_Q,Matrix(d'),n[2:end],k+1)')

    y = reshape(KPShiftSolve(TT,n,e,lambda,1.0),size(d))
    x = real(Matrix(v*kron_prod_times_matrix(Q,Matrix(y'),n[2:end],k+1)'))

    return x

end

function dlyap(a::Array{T,2}, b::Array{T,2}) where {T <: Real}

    n = size(a,1)
    x = zeros(n,n)
    j = n

    (s,u) = schur(a)
    b = u'b*u

    while j > 0
        j1 = j
        if j == 1
            block = 1
        elseif !isequal(s[j,j-1],0.0)
            block = 2
            j -= 1
        else
            block = 1
        end
        @views lhs = kron(s[j:j1,j:j1],s) - I # I = eye(block*n)
        @views rhs = vec(b[:,j:j1])
        if j1 < n
            @views rhs2 = s*(x[:,(j1+1):n]*s[j:j1,(j1+1):n]')
            rhs += vec(rhs2)
        end
        w = -lhs\rhs
        @views x[:,j] = w[1:n]
        if block == 2
            @views x[:,j1] = w[(n+1):block*n]
        end
        j -= 1
    end

    x = u*x*u'

    return x

end

function compute_variances(soln::FirstOrderSolutionStoch)

    hx = soln.hx
    k = soln.k
    gx = soln.gx
    sigma = soln.sigma

    var_states = dlyap(hx,k*sigma*k')
    var_jumps  = gx*var_states*gx'

  return var_states,var_jumps

end

function compute_chebyshev_integrals(eps_nodes,eps_weights,nodes,order,rho,sigma)

  # A simplified implementation of the integration described in Judd et al (2017, section 5.1);
  # it uses the mean integral across the nodes.
  # Using integrals[i] = sum(exp.(sqrt(2)*sigma*(i-1)*eps_nodes).*eps_weights)*pi^(-1/2)
  # leads to something much simplier and equally accurate, but is less recognizably appropriate
  # or the case where the shocks are AR(1) processes and ordinary polynomials are not being used.

  terms_num  = Array{Float64}(undef,length(eps_nodes))
  integrals  = Array{Float64}(undef,order+1,length(nodes))
  integrals2 = Array{Float64}(undef,order+1)
  for i = 1:(order+1)
    integrals2[i] = sum(exp.(sqrt(2)*sigma*(i-1)*eps_nodes).*eps_weights)*pi^(-1/2)
    for j = 1:length(nodes)
      terms_num     .= rho*nodes[j] .+ sqrt(2)*sigma*eps_nodes
      terms_den      = rho*nodes[j]
      terms_num     .= chebyshev_polynomial(i,terms_num)[:,i]
      terms_den      = chebyshev_polynomial(i,terms_den)[i]
      integrals[i,j] = sum((terms_num/terms_den).*eps_weights)*pi^(-1/2)
     end
  end

  nodetoosmall = abs.(nodes) .< sqrt(eps())
  if sum(nodetoosmall) > 0
    if length(nodes) == 1
      integrals[:,1] .= integrals2
    else
      for i = 1:length(nodes)
        if nodetoosmall[i] == 1
          if i == 1
            integrals[:,i] .= integrals[:,i+1]
        elseif i == length(nodes)
            integrals[:,i] .= integrals[:,i-1]
          else
            integrals[:,i] .= (integrals[:,i-1]+integrals[:,i+1])/2
          end
        end
      end
    end
  end

  return reshape(sum(integrals,dims=2)/length(nodes),order+1)

#  integrals = Array{Float64}(undef,order+1)
#  for i = 1:(order+1)
#    integrals[i] = sum(exp.(sqrt(2)*sigma*(i-1)*eps_nodes).*eps_weights)*pi^(-1/2)
#  end
#  return integrals

end

function compute_smolyak_integrals(eps_nodes,eps_weights,nx,order,grid,RHO,sigma)

  integrals = ones(nx,order+1)
  for j = 1:size(sigma,2)
      nodes = unique(grid[:,j])
      integrals[j,:] .= compute_chebyshev_integrals(eps_nodes,eps_weights,nodes,order,RHO[j,j],sigma[j,j])
  end

#  for j = 1:size(sigma,2)
#      for i = 1:(order+1)
#          integrals[j,i] = sum(exp.(sqrt(2)*sigma[j,j]*(i-1)*eps_nodes).*eps_weights)*pi^(-1/2)
#      end
#  end

  return integrals

end

function weight_scale_factors(eps_nodes,eps_weights,multi_index,nx,grid,RHO,sigma)

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i(unique_multi_index).-1

  # Here we construct the base integrals

  base_integrals = Array{Array{Float64,2}}(undef,length(unique_orders))
  for i = 1:length(unique_orders)
    base_integrals[i] = compute_smolyak_integrals(eps_nodes,eps_weights,nx,unique_orders[i],grid,RHO,sigma)
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_integrals = Array{Array{Float64,2}}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_integrals[i] = base_integrals[i][:,size(base_integrals[i-1],2)+1:end]
  end
  unique_base_integrals[1] = base_integrals[1]

  # Construct the first row of the interplation matrix

  new_integrals = unique_base_integrals[multi_index[1,1]][1,:]
  for i = 2:size(multi_index,2)
    new_integrals = kron(new_integrals,unique_base_integrals[multi_index[1,i]][i,:])
  end

  weight_scale_factor = copy(new_integrals)

  # Iterate over nodes, doing the above three steps at each iteration

  for j = 2:size(multi_index,1)

    new_integrals = unique_base_integrals[multi_index[j,1]][1,:]
    for i = 2:size(multi_index,2)
      new_integrals = kron(new_integrals,unique_base_integrals[multi_index[j,i]][i,:])
    end
    weight_scale_factor = [weight_scale_factor; new_integrals]

  end

  return weight_scale_factor

end

function scale_weights(weights,weight_scale_factor)

    scaled_weights = weights.*weight_scale_factor

    return scaled_weights

end

function compute_piecewise_linear_integrals(eps_nodes,eps_weights,sigma)

  integral = 1.0
  #integral = sum(exp.(sqrt(2)*sigma*eps_nodes).*eps_weights)*pi^(-1/2)

  return integral

end
