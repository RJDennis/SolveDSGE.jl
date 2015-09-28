function householder{T<:AbstractFloat}(a::Array{T,1})

  n = length(a)
  m = 1

  q = eye(n)
  r = copy(a)

  for i = 1:m

    x = r[i:n,i]
    e = zeros(n-i+1)
    e[1] = norm(x)
    u = x-e
    u = u/norm(u)
    qq = eye(n-i+1)-2*u*u'
    r[i:n,:] = qq*r[i:n,:]
    q[i:n,:] = qq*q[i:n,:]

  end

  return q',r

end

function  householder{T<:AbstractFloat}(a::Array{Complex{T},1})

  n = length(a)
  m = 1

  q = complex(eye(n))
  r = copy(a)

  for i = 1:m

    x = r[i:n,i]
    e = complex(zeros(n-i+1))
    e[1] = norm(x)
    u = x-e
    u = u/norm(u)
    qq = complex(eye(n-i+1))-(complex(1.0)+(x'u)/(u'x)).*u*u'
    r[i:n,:] = qq*r[i:n,:]
    q[i:n,:] = qq*q[i:n,:]

  end

  return q',r

end

function  householder{T<:AbstractFloat}(a::Array{T,2})

  (n,m) = size(a)

  q = eye(n)
  r = copy(a)

  for i = 1:m

    x = r[i:n,i]
    e = zeros(n-i+1)
    e[1] = norm(x)
    u = x-e
    u = u/norm(u)
    qq = eye(n-i+1)-2*u*u'
    r[i:n,:] = qq*r[i:n,:]
    q[i:n,:] = qq*q[i:n,:]

  end

  return q',r

end

function  householder{T<:AbstractFloat}(a::Array{Complex{T},2})

  (n,m) = size(a)

  q = complex(eye(n))
  r = copy(a)

  for i = 1:m

    x = r[i:n,i]
    e = complex(zeros(n-i+1))
    e[1] = norm(x)
    u = x-e
    u = u/norm(u)
    qq = complex(eye(n-i+1))-(complex(1.0)+(x'u)/(u'x)).*u*u'
    r[i:n,:] = qq*r[i:n,:]
    q[i:n,:] = qq*q[i:n,:]

  end

  return q',r

end
