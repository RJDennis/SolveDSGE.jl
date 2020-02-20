#=

Purpose: Use the Schur method to find the bounded solution of the
         discrete Lyapunov equation:

         X = B + A*X*A'.

Inputs:  A (n*n)
         B (n*n) = Symmetric positive semidefinite

Outputs: X (n*n) = Symmetric positive semidefinite solution.

Based on Octave code written by A.S. Hodel.

=#

function dlyap(a::Array{T, 2}, b::Array{T, 2}) where T <: AbstractFloat

  n = size(a,1)
  x = zeros(n,n)
  j = n

  (s,u) = schur(a)
  b = u'b*u

  while j > 0
    j1 = j
    if j == 1
      block = 1
    elseif s[j,j-1] > 2*eps(T)*norm(s)
      block = 2
      j -= 1
    else
      block = 1
    end
    lhs = kron(s[j:j1,j:j1],s) - I # I = eye(block*n)
    rhs = vec(b[:,j:j1])
    if j1 < n
      rhs2 = s*(x[:,(j1+1):n]*s[j:j1,(j1+1):n]')
      rhs += vec(rhs2)
    end
    w = -lhs\rhs
    x[:,j] = w[1:n]
    if block == 2
      x[:,j1] = w[(n+1):block*n]
    end
    j -= 1
  end

  x = u*x*u'

  return x

end
