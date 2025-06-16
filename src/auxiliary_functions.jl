########################## Auxiliary functions #############################

"""
Creates identity matrix.

Internal function; not exposed to users.
"""
function eye(tpe::Type,n::S) where {S<:Integer}

    if tpe <: Number

        return Matrix{tpe}(I,n,n)

    end

end

eye(n::Integer) = eye(Float64,n::Integer)

"""
Computes the matrix-trace as defined by Gomme and Klein (2011)

Internal function; not exposed to users.
"""
function tracem(x::Array{T,2}) where {T<:Real}

    (n1,n2) = size(x)

    if n1 == n2 # Matrix is square

        y = reshape([tr(x)],1,1)

        return y

    elseif n1 >= n2 # Matrix is skinny

        m = div(n1,n2)
        y = Array{T,2}(undef,m,1)  # We want this to be a 2-d array for subsequent matrix multiplication

        @inbounds for i = 1:m
            @views y[i,1] = tr(x[(n2*(i-1)+1):i*n2,1:n2])
        end

        return y

    else # Matrix is fat

        m = div(n2,n1)
        y = Array{T,2}(undef,1,m)  # We want this to be a 2-d array for subsequent matrix multiplication

        @inbounds for i = 1:m
            @views y[1,i] = tr(x[1:n1,(n1*(i-1)+1):i*n1])
        end

        return y

    end

end

"""
Efficiently computes (a ⊗ b) × v where a and b are matrices and v is a vector.

Internal function; not exposed to users.
"""
function kron_prod_times_vector(a::AbstractArray{T,2},b::AbstractArray{T,2},v::AbstractArray{T,1}) where {T<:Real}

    (n1, n2) = size(a)
    (n3, n4) = size(b)
    n5 = length(v)

    if n5 != n2*n4
        error("a, b, and v do not have the correct sizes.")
    end

    p = vec(b*reshape(v,n4,n2)*a')

    return p

end

"""
Efficiently computes v × (a ⊗ b) where a and b are matrices and v is a vector.

Internal function; not exposed to users.
"""
function vector_times_kron_prod(v::AbstractArray{T,1},a::AbstractArray{T,2},b::AbstractArray{T,2}) where {T<:Real}

    p = Matrix(kron_prod_times_vector(a',b',v')')

    return p

end

"""
Efficiently computes (a ⊗ b) × v where a and b and v are conformable matrices.

Internal function; not exposed to users.
"""
function kron_prod_times_matrix(a::AbstractArray{T,2},b::AbstractArray{T,2},v::AbstractArray{T,2}) where {T<:Real}

    (n1, n2) = size(a)
    (n3, n4) = size(b)
    (n5, n6) = size(v)

    if n5 != n2*n4
        error("a, b, and v do not have the correct sizes.")
    end

    p = zeros(n1*n3,n6)
    for i = 1:n6
        @views p[:,i] = kron_prod_times_vector(a,b,v[:,i])
    end

    return p

end

"""
Efficiently computes v × (a ⊗ b) where a and b and v are conformable matrices.

Internal function; not exposed to users.
"""
function matrix_times_kron_prod(v::AbstractArray{T,2},a::AbstractArray{T,2},b::AbstractArray{T,2}) where {T<:Real}

    p = Matrix(kron_prod_times_matrix(a',b',v')')

    return p

end

"""
Efficiently computes (a ⊗ b ⊗ c) × v where a, b, and c are matrices and v is a vector.

Internal function (not actually used); not exposed to users.
"""
function kron_prod_times_vector(a::AbstractArray{T,2},b::AbstractArray{T,2},c::AbstractArray{T,2},v::AbstractArray{T,1}) where {T<:Real}

    (n1, n2) = size(a)
    (n3, n4) = size(b)
    (n5, n6) = size(c)
    n7 = length(v)

    if n7 != n2* n4*n6
        error("a, b, c, and v do not have the correct sizes.")
    end

    v_tilda = reshape(v,n4*n6,n2)*a'
    p = vec(kron_prod_times_matrix(b,c,v_tilda))

    return p

end

"""
Efficiently computes v × (a ⊗ b ⊗ c) where a, b, and c # are matrices and v is a vector.

Internal function (not actually used); not exposed to users.
"""
function vector_times_kron_prod(v::AbstractArray{T,1},a::AbstractArray{T,2},b::AbstractArray{T,2},c::AbstractArray{T,2}) where {T<:Real}


    product = Matrix(kron_prod_times_vector(a',b',c',v')')

    return product

end

"""
Efficiently computes (a ⊗ b ⊗ c) × v where a, b, c and v are comformable matrices.

Internal function; not exposed to users.
"""
function kron_prod_times_matrix(a::AbstractArray{T,2},b::AbstractArray{T,2},c::AbstractArray{T,2},v::AbstractArray{T,2}) where {T<:Real}

    (n1, n2) = size(a)
    (n3, n4) = size(b)
    (n5, n6) = size(c)
    (n7, n8) = size(v)

    if n7 != n2*n4*n6
        error("a, b, c, and v do not have the correct sizes.")
    end

    p = zeros(n1*n3*n5,n8)
    for i = 1:n8
        @views v_tilda = reshape(v[:,i],n4*n6,n2)*a'
        p[:,i] .= vec(kron_prod_times_matrix(b,c,v_tilda))
    end

    return p

end

"""
Efficiently computes v × (a ⊗ b ⊗ c) where a, b, c and v are comformable matrices.

Internal function; not exposed to users.
"""
function matrix_times_kron_prod(v::AbstractArray{T,2},a::AbstractArray{T,2},b::AbstractArray{T,2},c::AbstractArray{T,2}) where {T<:Real}

    p = Matrix(kron_prod_times_matrix(a',b',c',v')')

    return p

end

"""
Efficiently computes (q1 ⊗ q2 ⊗ q3 ⊗ ... ⊗ qN) × v where the matrices in q and v are comformable.

Internal function; not exposed to users.
"""
function kron_prod_times_matrix(q::Array{Array{T,2},1},v::Array{T,2}) where {T<:Real}

    qnr = size.(q,1)
    qnc = size.(q,2)

    (vnr, vnc) = size(v)

    if vnr != prod(qnc)
        error("Matrices are not conformable.")
    end

    if length(q) == 2
        p = kron_prod_times_matrix(q[1],q[2],v)
        return p
    else
        p = zeros(prod(qnr),vnc)
        for i = 1:vnc
            @views v_tilda = reshape(v[:,i],prod(qnc[2:end]),qnc[1])*q[1]'
            p[:,i] .= vec(kron_prod_times_matrix(q[2:end],v_tilda))
        end
        return p
    end

end

"""
Efficiently computes v × (q1 ⊗ q2 ⊗ q3 ⊗ ... ⊗ qN) where the matrices in q and v are comformable.

Internal function; not exposed to users.
"""
function matrix_times_kron_prod(v::Array{T,2},q::Array{Array{T,2},1}) where {T<:Real}

    p = Matrix(kron_prod_times_matrix(reshape(Matrix.(q'),length(q)),Matrix(v'))')
    return p

end

"""
Uses the Hessenberg-Schur method to find the bounded solution of the
discrete Sylvester equation:

      X + A*X*B = C

Based on Golub, Nash, and Van Loan (1979).

Internal function; not exposed to users.
"""
function dsylvester(a::AbstractArray{T,2},b::AbstractArray{T,2},c::Union{AbstractArray{T,1},AbstractArray{T,2}}) where {T<:Real}

    n = size(a,1)
    m = size(b,1)
    x = zeros(size(c))

    (s, u) = schur(Matrix(b'))
    (v, t) = hessenberg(a)

    c = v'*c*u

    j = m
    while j > 0
        j1 = j
        if j == 1
            block_size = 1
        elseif isequal(s[j,j-1],0.0) == false
            block_size = 2
            j -= 1
        else
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

"""
Computes the matrix trace as defined by Binning (2013).  Used for # computing the second-order terms, hss, gss.

Internal function; not exposed to users.
"""
function trm(x::AbstractArray{T,2}) where {T<:Real}

    (n1, n2) = size(x)
    k = Int(round(sqrt(n2)))

    y = zeros(n1)
    for i = 1:k
        @views y += x[:,i+(i-1)*k]
    end

    return y

end

"""
Computes the matrix trace as defined by Binning (2013).  Used for computing the third-order terms, hssx, gssx.

Internal function; not exposed to users.
"""
function trm2(x::AbstractArray{T,2}) where {T<:Real}

    (n1, n2) = size(x)
    k = Int(round(n2^(1//3)))

    y = zeros(n1,k)
    for j = 1:k
        for i = 1:k
            @views y[:,j] += x[:,(j-1)+i+(i-1)*k^2]
        end
    end

    return y

end

"""
Computes the matrix trace.  Used for computing the fourth-order terms, hssxx, gssxx.

Internal function; not exposed to users.
"""
function trm3(x::AbstractArray{T,2}) where {T<:Real}

    (n1, n2) = size(x)
    k = Int(round(n2^(1//2)))

    y = zeros(n1,k)
    for j = 1:k
        for i = 1:k
            @views y[:,j] += x[:,j+(i-1)*k]
        end
    end

    return y

end

"""
Creates the combination matrix for a third-order perturbation as defined in Levintal (2017).

Internal function; not exposed to users
"""
function create_omega3(n::S) where {S<:Integer}

    # This function is a simplified version of the create_OMEGA function originally written in
    # Matlab by Oren Levintal for his paper "Fifth Order Perturbation Solution to DSGE Models"
    # published in the Journal of Economic Dynamics and Control, 2017.  Permission to translate
    # this function into Julia and release it within the SolveDSGE module was granted by Oren
    # Levintal on February 5, 2020.

    ind = [1:n^3;]
    M = reshape(ind,1,n,n,n)
    Ix = eye(S,n^3)
    omega1 = (reshape(Ix[:,PermutedDimsArray(M,[1, 4, 2, 3])],n^3,n^3)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 4, 3])],n^3,n^3)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 4])],n^3,n^3))

    return omega1

end

"""
Creates the combination matrices for a fourth-order perturbation as defined in Levintal (2017).

Internal function; not exposed to users.
"""
function create_omega4(n::S) where {S<:Integer}

    # This function is a simplified version of the create_OMEGA function originally written in
    # Matlab by Oren Levintal for his paper "Fifth Order Perturbation Solution to DSGE Models"
    # published in the Journal of Economic Dynamics and Control, 2017.  Permission to translate
    # this function into Julia and release it within the SolveDSGE module was granted by Oren
    # Levintal on February 5, 2020.

    ind = [1:n^4;]
    M = reshape(ind,1,n,n,n,n)
    Ix = eye(S,n^4)
    omega2 = (reshape(Ix[:,PermutedDimsArray(M,[1, 5, 4, 2, 3])],n^4,n^4)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 5, 2, 4, 3])],n^4,n^4)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 5, 4, 3])],n^4,n^4)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 5, 2, 3, 4])],n^4,n^4)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 5, 3, 4])],n^4,n^4)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 5, 4])],n^4,n^4))

    omega3 = (reshape(Ix[:,PermutedDimsArray(M,[1, 5, 2, 3, 4])],n^4,n^4)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 5, 3, 4])],n^4,n^4)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 5, 4])],n^4,n^4)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 4, 5])],n^4,n^4))

    omega4 = (reshape(Ix[:,PermutedDimsArray(M,[1, 4, 2, 3, 5])],n^4,n^4)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 4, 3, 5])],n^4,n^4)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 4, 5])],n^4,n^4))

    return omega2, omega3, omega4

end

"""
Creates the combination matrices for a fifth-order perturbation as defined in Levintal (2017).

Internal function (not actually used); not exposed to users.
"""
function create_omega5(n::S) where {S<:Integer}

    # This function is a simplified version of the create_OMEGA function originally written in
    # Matlab by Oren Levintal for his paper "Fifth Order Perturbation Solution to DSGE Models"
    # published in the Journal of Economic Dynamics and Control, 2017.  Permission to translate
    # this function into Julia and release it within the SolveDSGE module was granted by Oren
    # Levintal on February 5, 2020.

    ind = [1:n^5;]
    M = reshape(ind,1,n,n,n,n,n)
    Ix = eye(S,n^5)
    omega5 = (reshape(Ix[:,PermutedDimsArray(M,[1, 6, 5, 4, 2, 3])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 6, 5, 2, 4, 3])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 6, 2, 5, 4, 3])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 6, 5, 4, 3])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 6, 5, 2, 3, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 6, 2, 5, 3, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 6, 5, 3, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 6, 2, 3, 5, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 6, 3, 5, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 6, 5, 4])],n^5,n^5))

    omega6 = (reshape(Ix[:,PermutedDimsArray(M,[1, 6, 5, 2, 3, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 6, 2, 5, 3, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 6, 5, 3, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 6, 2, 3, 5, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 6, 3, 5, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 6, 5, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 6, 2, 3, 4, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 6, 3, 4, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 6, 5, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 4, 6, 5])],n^5,n^5))

    omega7 = (reshape(Ix[:,PermutedDimsArray(M,[1, 6, 4, 2, 3, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 6, 2, 4, 3, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 6, 4, 3, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 6, 2, 3, 4, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 6, 3, 4, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 6, 4, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 4, 6, 2, 3, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 4, 2, 6, 3, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 4, 6, 3, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 4, 2, 3, 6, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 4, 3, 6, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 4, 6, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 4, 2, 3, 6, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 4, 3, 5, 6])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 4, 5, 6])],n^5,n^5))

    omega8 = (reshape(Ix[:,PermutedDimsArray(M,[1, 6, 2, 3, 4, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 6, 3, 4, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 6, 4, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 4, 6, 5])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 4, 5, 6])],n^5,n^5))

    omega9 = (reshape(Ix[:,PermutedDimsArray(M,[1, 5, 2, 3, 4, 6])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 5, 3, 4, 6])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 5, 4, 6])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 4, 5, 6])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 5, 2, 3, 6, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 5, 3, 6, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 3, 5, 6, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 2, 5, 6, 3, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 5, 2, 6, 3, 4])],n^5,n^5)
            + reshape(Ix[:,PermutedDimsArray(M,[1, 5, 6, 2, 3, 4])],n^5,n^5))

    return omega5, omega6, omega7, omega8, omega9

end

"""
Computes y = (A[p]⊗A[p-1]⊗...⊗A[1] )*x where the elements of A and matrices and x is a conformable vector.

Internal function; not exposed to users.
"""
function kron_prod_times_vector(A::Union{Array{Array{T,2},1},Array{Array{Complex{T},2},1}},x::Union{Array{T,1},Array{Complex{T},1}},n::Array{S,1},p::S) where {T<:Real,S<:Integer}

    N = prod(n[1:p])
    z = copy(x)
    for i = 1:p
        z = transpose(A[i]*reshape(z,n[i],div(N,n[i])))
    end
    y = reshape(Matrix(z),N)

    return y

end

"""
Computes y = (A[p]⊗A[p-1]⊗...⊗A[1] )*x where the elements of A and x are conformable matrices.

Internal function; not exposed to users.
"""
function kron_prod_times_matrix(A::Union{Array{Array{T,2},1},Array{Array{Complex{T},2},1}},x::Union{Array{T,2},Array{Complex{T},2}},n::Array{S,1},p::S) where {T<:Real,S<:Integer}

    N = prod(n[1:p])
    y = Array{Complex{T}}(undef,N,size(x,2))
    for j in axes(x,2)
        @views z = x[:,j]
        for i = 1:p
            z = transpose(A[i]*reshape(z,n[i],div(N,n[i])))
        end
        y[:,j] .= reshape(Matrix(z),N)
    end

    return y

end

function KPShiftSolve(TT::Union{Array{Array{T,2},1},Array{Array{Complex{T},2},1}},n::Array{S,1},c::Union{Array{T,1},Array{Complex{T},1}},lambda::T,alpha::Union{T,Complex{T}}) where {T<:Real,S<:Integer}

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

"""
Uses a recursive Schur method to find the bounded solution of the Sylvester
equation:

AX + B*X*(Kron^(k)C) = D

Based on Martin and Van Loan (2006).  This is a simplified implementation
of their algorithm, but it captures most of the gains over Golub, Nash,
and van Loan (1979).

Internal function; not exposed to users.
"""
function martin_van_loan(a::Array{T,2},b::Array{T,2},c::Array{T,2},d::Array{T,2},k::S) where {T<:Real,S<:Integer}

    a = copy(a)
    b = a\copy(b)
    c = copy(c)
    d = a\copy(d)

    (v, s) = hessenberg(b)       # v*s*v' = b
    (t, q) = schur(complex(c'))  # q*t*q' = c'

    v = Matrix(v)

    p = k + 2
    TT = Array{typeof(t)}(undef,p)
    TT[1] = s
    for i = 2:p
        TT[i] = conj(t)
    end

    Q = fill(q,k + 1)
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

"""

Solves discrete Lyapunov equations.

Internal function; not exposed to users.
"""
function dlyap(a::Array{T,2},b::Array{T,2}) where {T<:Real}

    n = size(a,1)
    x = zeros(n,n)
    j = n

    (s, u) = schur(a)
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

"""
Converts from a Cartesian Index to subscripts.

Internal function; not exposed to users.
"""
function ind2sub(i::S,dims::Tuple{S,Vararg{S}}) where {S<:Integer}

    if i < 1 || i > prod(dims)
        error("Index is out of bounds.")
    end

    subs = CartesianIndices(dims)[i]

    return subs

end

"""
Computes the unconditional variance of a first-order VAR.

Internal function; not exposed to users.
"""
function compute_variances(soln::FirstOrderSolutionStoch)

    hx = soln.hx
    k = soln.k
    gx = soln.gx
    sigma = soln.sigma

    var_states = dlyap(hx,k*sigma*k')
    var_jumps = gx*var_states*gx'

    return var_states, var_jumps

end

"""
Computes the integrals needed for precomputed quadrature for Chebyshev polynomials 
for the case where the shocks are AR(1) and the innovations are independent.
    
Internal function; not exposed to users.
"""
function _compute_chebyshev_integrals(eps_nodes::Array{T,1},eps_weights::Array{T,1},nodes::Array{T,1},order::S,rho::T,sigma::T) where {T<:AbstractFloat,S<:Integer}

    # Case where shocks are AR(1) and innovations are independent

    if rho == 0.0
        rho = eps()
    end

    selected_nodes = copy(nodes)
    if length(nodes) > 2
        selected_number = ceil(Int,length(nodes)^(1/2))
        if isodd(length(nodes)) && iseven(selected_number)
            selected_number -= 1
        elseif iseven(length(nodes)) && isodd(selected_number)
            selected_number -= 1
        end
        start_index = Int((length(nodes) - selected_number)/2)
        selected_nodes = nodes[start_index+1:start_index+selected_number]
    end

    terms_num = Array{T}(undef,length(eps_nodes))
    integrals = Array{T,2}(undef,order+1,length(selected_nodes))
    integrals2 = Array{T}(undef,order+1)
    for i = 1:(order+1)
        integrals2[i] = sum(exp.(sqrt(2)*sigma*(i-1)*eps_nodes).*eps_weights)*pi^(-1/2)
        for j in eachindex(selected_nodes)
            terms_num .= rho*nodes[j] .+ sqrt(2)*sigma*eps_nodes
            terms_den = rho*nodes[j]
            terms_num .= chebyshev_polynomial(i,terms_num)[:,i]
            terms_den = chebyshev_polynomial(i,terms_den)[i]
            integrals[i,j] = sum((terms_num/terms_den).*eps_weights)*pi^(-1/2)
        end
    end

    nodetoosmall = abs.(selected_nodes) .< sqrt(eps())
    if sum(nodetoosmall) > 0
        if length(selected_nodes) == 1
            integrals[:,1] .= integrals2
        else
            for i in eachindex(selected_nodes)
                if nodetoosmall[i] == 1
                    if i == 1
                        integrals[:,i] .= integrals[:,i+1]
                    elseif i == length(selected_nodes)
                        integrals[:,i] .= integrals[:,i-1]
                    else
                        integrals[:,i] .= (integrals[:,i-1] + integrals[:,i+1])/2
                    end
                end
            end
        end
    end

    return reshape(sum(integrals,dims = 2)/length(selected_nodes),order+1)

end

"""
Computes the integrals needed for precomputed quadrature for Chebyshev polynomials 
for the case where the shocks are AR(1) and the innovations are correlated.

Internal function; not exposed to users.
"""
function _compute_chebyshev_integrals(eps_nodes::Array{T,1},eps_weights::Array{T,1},nodes::Array{Array{T,1},1},order::Union{S,Array{S,1}},Ρ::Array{T,2},k::Array{T,2}) where {T<:AbstractFloat,S<:Integer}

    # Case where shocks are AR(1) and innovations are correlated

    n = size(k,2)
    for i = 1:n
        if Ρ[i,i] == 0.0
            Ρ[i,i] = eps()
        end
    end

    selected_nodes = similar(nodes)
    for i in eachindex(nodes)
        if length(nodes[i]) > 2
            selected_number = ceil(Int,sqrt(length(nodes[i])))
            if isodd(length(nodes[i])) && iseven(selected_number)
                selected_number -= 1
            elseif iseven(length(nodes[i])) && isodd(selected_number)
                selected_number -= 1
            end
            start_index = Int((length(nodes[i]) - selected_number)/2)
            selected_nodes[i] = nodes[i][start_index+1:start_index+selected_number]
        end
    end

    if typeof(order) == S
        ord = Tuple(fill(order,N))
    else
        ord = Tuple(order)
    end

    terms_num = Array{T}(undef,length(eps_nodes))
    integrals = Array{T}(undef,(ord.+1...,length.(selected_nodes)...))
    integrals2 = Array{T}(undef,ord.+1)

    order_prod = prod(ord.+1)
    nodes_prod = prod(length.(selected_nodes))
    eps_prod = length(eps_nodes)^n

    for i = 1:order_prod
        ii = ind2sub(i,ord.+1)
        int2 = 0.0
        for j = 1:eps_prod
            jj = ind2sub(j,Tuple(fill(length(eps_nodes),n)))
            eps_w = eps_weights[jj[1]]
            eps_node = eps_nodes[jj[1]]
            for k = 2:n
                eps_w *= eps_weights[jj[k]]
                eps_node = [eps_node; eps_nodes[jj[k]]]
            end
            int2 += exp(2^(n/2)*sqrt(det(k*k'))*collect(Tuple(ii).-1)'*(k*k')*eps_node)*eps_w*pi^(-n/2)
        end
        integrals2[Tuple(ii)...] = int2 # exp(collect(ii.-1)'*(k*k')*collect(ii.-1)/2) # This is the analytic expression
    end

    for i = 1:order_prod
        ii = ind2sub(i,ord.+1)
        for m = 1:nodes_prod
            mm = ind2sub(m,Tuple(length.(selected_nodes)))
            int2 = 0.0
            node = selected_nodes[1][mm[1]]
            for k = 2:n
                node = [node; selected_nodes[k][mm[k]]]
            end
            for j = 1:eps_prod
                jj = ind2sub(j,Tuple(fill(length(eps_nodes),n)))
                eps_w = eps_weights[jj[1]]
                eps_node = eps_nodes[jj[1]]
                for k = 2:n
                    eps_w *= eps_weights[jj[k]]
                    eps_node = [eps_node; eps_nodes[jj[k]]]
                end
                terms_num = Ρ*node + sqrt(2)*k*eps_node
                terms_den = Ρ*node
                num = chebyshev_polynomial(ii[1],terms_num[1])[ii[1]]
                den = chebyshev_polynomial(ii[1],terms_den[1])[ii[1]]
                for k = 2:n
                    num = [num; chebyshev_polynomial(ii[k],terms_num[k])[ii[k]]]
                    den = [den; chebyshev_polynomial(ii[k],terms_den[k])[ii[k]]]
                end
                int2 += (prod(num)/prod(den))*eps_w*pi^(-n/2)
            end
            integrals[Tuple(ii)...,Tuple(mm)...] = int2
        end
    end

    for i in eachindex(integrals)
        ii = ind2sub(i,Tuple(size(integrals)))
        if abs(integrals[i]) > 2.0 || isnan(integrals[i])
            if length(integrals) == 1
                integrals[i] = integrals2[i]
            else
                integrals[i] = integrals2[CartesianIndex(Tuple(ii)[1:n])]
            end
        end
    end

    for i = n:-1:1
        integrals = sum(integrals,dims = (n+i))/length(selected_nodes[i])
    end

    return reshape(integrals,ord.+1)

end

"""
Computes the integrals needed for precomputed quadrature with Chebyshev plynomials.

Internal function; not exposed to users.
"""
function compute_chebyshev_integrals(eps_nodes::Array{T,1},eps_weights::Array{T,1},nodes::Array{Array{T,1},1},order::Union{S,Array{S,1}},Ρ::Array{T,2},k::Array{T,2}) where {T<:AbstractFloat,S<:Integer}

    if !isdiag(Ρ .> sqrt(eps()))
        error("The autoregression matrix for the shocks must be diagonal.")
    end

    ns = size(k,2)

    if typeof(order) == S
        ord = fill(order,ns)
    else
        ord = order[1:ns]
    end

    if !isdiag(abs.(k) .> sqrt(eps())) # Correlated innovations
        integrals = _compute_chebyshev_integrals(eps_nodes,eps_weights,nodes[1:ns],ord,Ρ,k)
        return integrals
    else # Uncorrelated innovations
        integrals = Array{Array{T,1},1}(undef,ns)
        for i = 1:ns
            integrals[i] = _compute_chebyshev_integrals(eps_nodes,eps_weights,nodes[i],ord[i],Ρ[i,i],k[i,i])
        end
        return integrals
    end

end

"""
Uses the precomputed integrals to scale the Chebyshev weights.

Internal function; not exposed to users.
"""
function scale_chebyshev_weights!(weights::Array{Array{T,N},1},scaled_weights::Array{Array{T,N},1},integrals::Array{Array{T,1},1},j_approx::Union{S,Array{S,1}},ns::S) where {T<:AbstractFloat,N,S<:Integer}

    for i in eachindex(j_approx)
        for j = 1:ns
            index = [1:ndims(weights[i]);]
            index[1],index[j] = index[j],index[1]
            scaled_weights[i] .= permutedims(integrals[j].*permutedims(weights[i],index),index)
        end
    end

end

function scale_chebyshev_weights!(weights::Array{Array{T,N},1},scaled_weights::Array{Array{T,N},1},integrals::Array{T,N2},j_approx::Union{S,Array{S,1}},ns::S) where {T<:AbstractFloat,N,N2,S<:Integer}

    for i in eachindex(j_approx)
        for j = 1:ns
            scaled_weights[i] .= integrals.*weights[i]
        end
    end

end

"""
Computes the integral factors for Smolyak and the hyperbolic cross approximation for the 
case where the shocks are AR(1) and the innovations are independent.

Internal function; not exposed to users.
"""
function _compute_sparse_integrals(eps_nodes::Array{T,1},eps_weights::Array{T,1},nodes::Array{T,1},order::S,rho::T,sigma::T) where {T<:AbstractFloat,S<:Integer}

    # Case where shocks are AR(1) and innovations are independent

    if rho == 0.0
        rho = eps()
    end

    terms_num = Array{T}(undef,length(eps_nodes))
    integrals = Array{T,2}(undef,order+1,length(nodes))
    integrals2 = Array{T}(undef,order+1)
    for i = 1:(order+1)
        integrals2[i] = sum(exp.(sqrt(2)*sigma*(i-1)*eps_nodes).*eps_weights)*pi^(-1/2)
        for j in eachindex(nodes)
            terms_num .= rho*nodes[j] .+ sqrt(2)*sigma*eps_nodes
            terms_den = rho*nodes[j]
            terms_num .= chebyshev_polynomial(i,terms_num)[:,i]
            terms_den = chebyshev_polynomial(i,terms_den)[i]
            integrals[i,j] = sum((terms_num/terms_den).*eps_weights)*pi^(-1/2)
        end
    end

    nodetoosmall = abs.(nodes) .< sqrt(eps())
    if sum(nodetoosmall) > 0
        if length(nodes) == 1
            integrals[:,1] .= integrals2
        else
            for i in eachindex(nodes)
                if nodetoosmall[i] == 1
                    if i == 1
                        integrals[:,i] .= integrals[:,i+1]
                    elseif i == length(nodes)
                        integrals[:,i] .= integrals[:,i-1]
                    else
                        integrals[:,i] .= (integrals[:,i-1] + integrals[:,i+1])/2
                    end
                end
            end
        end
    end

    return reshape(sum(integrals,dims = 2)/length(nodes),order+1)

end

"""
Computes the integral factors needed for Smolyak and hyperbolic cross approximation.

Internal function; not exposed to users.
"""
function compute_sparse_integrals(eps_nodes,eps_weights,nx,order,grid,RHO,k)

    integrals = ones(nx,order+1)
    for j in axes(k,2)
        nodes = unique(grid[:,j])
        integrals[j,:] .= _compute_sparse_integrals(eps_nodes,eps_weights,nodes,order,RHO[j,j],k[j,j])
    end

    return integrals

end

"""
Uses the precomputed integrals to construct the scale factors for Smolyak polynomials.

Internal function; not exposed to users.
"""
function smolyak_weight_scale_factors(eps_nodes,eps_weights,multi_index,nx,grid,RHO,sigma)

    unique_multi_index = sort(unique(multi_index))
    unique_orders = SmolyakApprox.m_i(unique_multi_index) .- 1

    # Here we construct the base integrals

    base_integrals = Array{Array{Float64,2}}(undef,length(unique_orders))
    for i in eachindex(unique_orders)
        base_integrals[i] = compute_sparse_integrals(eps_nodes,eps_weights,nx,unique_orders[i],grid,RHO,sigma)
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

"""
Uses the precomputed integrals to construct the scale factors for hyperbolic cross polynomials.

Internal function; not exposed to users.
"""
function hcross_weight_scale_factors(eps_nodes,eps_weights,multi_index,nx,grid,RHO,sigma)

    T = eltype(eps_nodes)

    n = 2*maximum(multi_index,dims = 1) .+ 1
    (N, d) = HyperbolicCrossApprox.determine_grid_size(multi_index)

    base_integrals = Array{Array{T,2},1}(undef,d)
    unique_base_integrals = Array{Array{T,1},2}(undef,size(multi_index))

    # Construct the base polynomials

    for i = 1:d
        base_integrals[i] = compute_sparse_integrals(eps_nodes,eps_weights,nx,n[i]-1,grid,RHO,sigma)
    end

    # Compute the unique polynomial terms from the base polynomials

    for i in axes(multi_index,1)
        for j = 1:d
            if multi_index[i,j] == 0
                unique_base_integrals[i,j] = [base_integrals[j][1]]
            else
                unique_base_integrals[i,j] = [base_integrals[j][2*multi_index[i,j]],base_integrals[j][2*multi_index[i,j]+1]]
            end
        end
    end

    # Construct the first row of the interplation matrix

    weight_scale_factor = Array{T,1}(undef,N)
    l = 1
    @inbounds for j in axes(multi_index,1)
        new_integrals = unique_base_integrals[j,1]
        for i = 2:d
            new_integrals = kron(new_integrals,unique_base_integrals[j,i])
        end
        m = length(new_integrals)
        weight_scale_factor[l:l+m-1] = new_integrals
        l += m
    end

    return weight_scale_factor

end

"""
Uses the precomputed scale factors to scale the weight in either Smolyak or hyperbolic cross polynomials.

Internal function; not exposed to users.
"""
function scale_sparse_weights(weights,weight_scale_factor)

    scaled_weights = weights.*weight_scale_factor

end

"""
Computes the integrals needed for precomputed quadrature using piecewise linear approximation.

Internal function; not exposed to users.
"""
function compute_piecewise_linear_integrals(eps_nodes,eps_weights,sigma)

    integral = 1.0
    #integral = sum(exp.(sqrt(2)*sigma*eps_nodes).*eps_weights)*pi^(-1/2)

    return integral

end

"""
Computes the maximum absolute distance between elements in two multi-dimensional arrays.

Internal function; not exposed to users.
"""
function mylen(len::T,a::Array{T,N},b::Array{T,N}) where {T<:AbstractFloat,N}

    for i in eachindex(a)
        len = max(len,abs(a[i]-b[i]))
    end

    return len

end