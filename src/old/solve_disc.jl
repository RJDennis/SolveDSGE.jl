function solve_disc(model::State_Space_Form{T}, obj::State_Space_Objective{T}, tol::T, maxiters::S) where {T <: AbstractFloat, S <: Int}

  nx    = copy(model.nx)
  ny    = copy(model.ny)
  beta  = copy(obj.beta)
  sigma = copy(model.sigma)

  a = copy(model.a)
  b = copy(model.b)
  c = copy(model.c)

  q = copy(obj.q)
  u = copy(obj.u)
  r = copy(obj.r)

  n = nx+ny

  a11 = a[1:nx,1:nx]
  a12 = a[1:nx,(nx+1):n]
  a21 = a[(nx+1):n,1:nx]
  a22 = a[(nx+1):n,(nx+1):n]

  b1 = b[1:nx,:]
  b2 = b[(nx+1):n,:]

  q11 = q[1:nx,1:nx]
  q12 = q[1:nx,(nx+1):n]
  q21 = q[(nx+1):n,1:nx]
  q22 = q[(nx+1):n,(nx+1):n]

  u1 = u[1:nx,:]
  u2 = u[(nx+1):n,:]

  j = (I-a22)\a21   # I = eye(ny)
  k = (I-a22)\b2    # I = eye(ny)
  f = (k'q22*k+k'u2+u2'k+r)\(k'q21+k'q22*j+u1'+u2'j)
  h = j-k*f
  v = q11+2*q12*h+h'q22*h-2*u1*f-2*h'u2*f+f'r*f

  len = Inf
  iters = 0
  while len > tol && iters < maxiters

    j    = (a22-h*a12)\(h*a11-a21)
    k    = (a22-h*a12)\(h*b1-b2)
    abar = a11+a12*j
    bbar = b1+a12*k
    qbar = q11+j'q21+q12*j+j'q22*j
    ubar = q12*k+j'q22*k+u1+j'u2
    rbar = r+k'q22*k+u2'k+k'u2
    fn   = (rbar+beta*bbar'v*bbar)\(ubar'+beta*bbar'v*abar)
    vn   = beta*(abar-bbar*fn)'v*(abar-bbar*fn)+qbar-fn'ubar'-ubar*fn+fn'rbar*fn
    hn   = j-k*fn

    len = maximum(maximum(abs,vn-v), maximum(abs,vn-v), maximum(abs,fn-f))
    h = hn
    v = vn
    f = fn

    iters += 1

  end

  if iters <= maxiters
    retcode = true
  else
    retcode = false
  end

  f = -f
  p = a11+a12*h+b1*f
  h = [eye(nx); h]
  k = c[1:nx,:]

  soln = State_Space_Soln(p,k,h,f,v,sigma,retcode)

  return soln

end

function solve_disc(model::Generalized_State_Space_Form{T}, obj::State_Space_Objective{T}, tol::T, maxiters::S) where {T <: AbstractFloat, S <: Int}

  nx    = copy(model.nx)
  ny    = copy(model.ny)
  beta  = copy(obj.beta)
  sigma = copy(model.sigma)

  a0 = copy(model.a0)
  a  = copy(model.a)
  b  = copy(model.b)
  c  = copy(model.c)

  q = copy(obj.q)
  u = copy(obj.u)
  r = copy(obj.r)

  n = nx+ny

  a11 = a[1:nx,1:nx]
  a12 = a[1:nx,(nx+1):n]
  a21 = a[(nx+1):n,1:nx]
  a22 = a[(nx+1):n,(nx+1):n]

  b1 = b[1:nx,:]
  b2 = b[(nx+1):n,:]

  q11 = q[1:nx,1:nx]
  q12 = q[1:nx,(nx+1):n]
  q21 = q[(nx+1):n,1:nx]
  q22 = q[(nx+1):n,(nx+1):n]

  u1 = u[1:nx,:]
  u2 = u[(nx+1):n,:]

  j = (a0-a22)\a21
  k = (a0-a22)\b2
  f = (k'q22*k+k'u2+u2'k+r)\(k'q21+k'q22*j+u1'+u2'j)
  h = j-k*f
  v = q11+2*q12*h+h'q22*h-2*u1*f-2*h'u2*f+f'r*f

  len = Inf
  iters = 0
  while len > tol && iters < maxiters

    j    = (a22-a0*h*a12)\(a0*h*a11-a21)
    k    = (a22-a0*h*a12)\(a0*h*b1-b2)
    abar = a11+a12*j
    bbar = b1+a12*k
    qbar = q11+j'q21+q12*j+j'q22*j
    ubar = q12*k+j'q22*k+u1+j'u2
    rbar = r+k'q22*k+u2'k+k'u2
    fn   = (rbar+beta*bbar'v*bbar)\(ubar'+beta*bbar'v*abar)
    vn   = beta*(abar-bbar*fn)'v*(abar-bbar*fn)+qbar-fn'ubar'-ubar*fn+fn'rbar*fn
    hn   = j-k*fn

    len = maximum(maximum(abs,vn-v), maximum(abs,vn-v), maximum(abs,fn-f))
    h = hn
    v = vn
    f = fn

    iters += 1

  end

  if iters <= maxiters
    retcode = true
  else
    retcode = false
  end

  f = -f
  p = a11+a12*h+b1*f
  h = [eye(nx); h]
  k = c[1:nx,:]

  soln = State_Space_Soln(p,k,h,f,v,sigma,retcode)

  return soln

end

function solve_disc(model::Structural_Form{T}, obj::Structural_Objective{T}, tol::T, maxiters::S) where {T <: AbstractFloat, S <: Int}

  a0    = copy(model.a0)
  a1    = copy(model.a1)
  a2    = copy(model.a2)
  a3    = copy(model.a3)
  a5    = copy(model.a5)
  sigma = copy(model.sigma)

  q = copy(obj.q)
  r = copy(obj.r)

  beta = copy(obj.beta)

  h1 = a0\a1
  h2 = a0\a5
  f1 = zeros(size(a3'))
  f2 = zeros(size(a3,2),size(a5,2))
  v = q

  a3_bar = Array(T,size(a3))
  a5_bar = Array(T,size(a5))

  len = Inf
  iters = 0
  while len > tol && iters < maxiters

    a1_bar = (a0-a2*h1)\a1
    a3_bar = (a0-a2*h1)\a3
    a5_bar = (a0-a2*h1)\a5

    vn = q+beta*(a1_bar+a3_bar*f1)'v*(a1_bar+a3_bar*f1)

    f1n = -(r+a3_bar'vn*a3_bar)\(a3_bar'vn*a1_bar)
    h1n = a1_bar+a3_bar*f1

    len = maximum(maximum(abs,vn-v), maximum(abs,h1n-h1), maximum(abs,f1n-f1))
    h1 = h1n
    f1 = f1n
    v  = vn

    iters += 1

  end

  if iters <= maxiters
    retcode = true
  else
    retcode = false
  end

  f2 = -(r+a3_bar'v*a3_bar)\(a3_bar'v*a5_bar)
  h2 = a5_bar+a3_bar*f2

  h = [[h1 zeros(size(a3))];
       [f1 zeros(size(r))]]

  g = [h2; f2]

  soln = Structural_Soln(h,g,v,sigma,retcode)

  return soln

end

function solve_disc(model::Generalized_Structural_Form{T}, obj::Structural_Objective{T}, tol::T, maxiters::S) where {T <: AbstractFloat, S <: Int}

  a0    = copy(model.a0)
  a1    = copy(model.a1)
  a2    = copy(model.a2)
  a3    = copy(model.a3)
  a4    = copy(model.a4)
  a5    = copy(model.a5)
  sigma = copy(model.sigma)

  q = copy(obj.q)
  r = copy(obj.r)

  beta = copy(obj.beta)

  h1 = a0\a1
  f1 = zeros(size(a3'))
  di = Array(T,size(a0))
  v  = q

  len = Inf
  iters = 0
  while len > tol && iters < maxiters

    (v,retcode_doubling) = doubling(beta*h1',h1,q+beta*f1'r*f1,tol,maxiters)
    di  = inv(a0-a2*h1-a4*f1)
    f1n = -(r+a3'di'v*di*a3)\(a3'di'v*di*a1)
    h1  = di*(a1+a3*f1n)

    len = maximum(abs,f1n-f1)
    f1 = f1n

    iters += 1

  end

  if iters <= maxiters
    retcode = true
  else
    retcode = false
  end

  f2 = -(r+a3'di'v*di*a3)\(a3'di'v*di*a5)
  h2 = (a0-a2*h1-a4*f1)\(a5+a3*f2)

  h = [h1 zeros(size(a3));
       f1 zeros(size(r))]

  g = [h2; f2]

  soln = Structural_Soln(h,g,v,sigma,retcode)

end
