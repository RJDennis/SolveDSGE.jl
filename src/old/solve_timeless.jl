function solve_timeless(model::State_Space_Form{T}, obj::State_Space_Objective{T}, tol::T, maxiters::S) where {T <: AbstractFloat, S <: Int}

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
  np = size(b,2)
  ns = size(c,2)

  v = q+(a+a')/n
  f = zeros(np,n)

  len = Inf
  iters = 0
  while len > tol && iters < maxiters

    vn = q-u*f-f'u'+f'r*f+beta*(a-b*f)'v*(a-b*f)
    fn = (r+beta*b'vn*b)\(u'+beta*b'vn*a)

    len = maximum(maximum(abs,vn-v), maximum(abs,fn-f))
    f = fn
    v = vn

    iters += 1

  end

  if iters <= maxiters
    retcode = true
  else
    retcode = false
  end

  t = eye(n)
  t[(nx+1):n,1:nx]     = v[(nx+1):n,1:nx]
  t[(nx+1):n,(nx+1):n] = v[(nx+1):n,(nx+1):n]

  p = t*(a-b*f)/t
  k = t*[c[1:nx,:]; -(v[(nx+1):n,(nx+1):n])\v[(nx+1):n,1:nx]*c[1:nx,:]]
  h = inv(t)
  v = h'v*h
  f = -f*h

  # Now re-express the solution by eliminating the Lagrange multipliers

  p11 = p[1:nx,1:nx]
  p12 = p[1:nx,(nx+1):n]
  p21 = p[(nx+1):n,1:nx]
  p22 = p[(nx+1):n,(nx+1):n]

  k11 = k[1:nx,:]

  h21 = h[(nx+1):n,1:nx]
  h22 = h[(nx+1):n,(nx+1):n]

  f1 = f[1:np,1:nx]
  f2 = f[1:np,(nx+1):n]

  l = [h22; f2]
  g = [h21; f1]

  invl = pinv(l)

  p0                             = eye(2*nx+ny+np)
  p0[1:nx,(2*nx+1):(2*nx+ny+np)] = -p12*invl

  p1                                              = zeros(2*nx+ny+np,2*nx+ny+np)
  p1[1:nx,1:nx]                                   = p11-p12*invl*g
  p1[(nx+1):(2*nx),1:nx]                          = eye(nx)
  p1[(2*nx+1):(2*nx+ny+np),1:nx]                  = g
  p1[(2*nx+1):(2*nx+ny+np),(nx+1):(2*nx)]         = l*(p21-p22*invl*g)
  p1[(2*nx+1):(2*nx+ny+np),(2*nx+1):(2*nx+ny+np)] = l*p22*invl

  p = p0\p1

  h = p[(nx+1):(2*nx+ny),:]
  f = p[(2*nx+ny+1):(2*nx+ny+np),:]
  k = [k11; zeros((nx+ny+np),ns)]

  soln = State_Space_Soln(p,k,h,f,v,sigma,retcode)

  return soln

end

function solve_timeless(model::Generalized_State_Space_Form{T}, obj::State_Space_Objective{T}, tol::T, maxiters::S) where {T <: AbstractFloat, S <: Int}

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

  n  = nx+ny
  np = size(b,2)
  ns = size(c,2)

  q11 = q[1:nx,1:nx]
  q12 = q[1:nx,(nx+1):n]
  q21 = q[(nx+1):n,1:nx]
  q22 = q[(nx+1):n,(nx+1):n]

  u1 = u[1:nx,:]
  u2 = u[(nx+1):n,:]

  a11 = a[1:nx,1:nx]
  a12 = a[1:nx,(nx+1):n]
  a21 = a[(nx+1):n,1:nx]
  a22 = a[(nx+1):n,(nx+1):n]

  b1 = b[1:nx,:]
  b2 = b[(nx+1):n,:]

  a_tilda = [[a11 zeros(nx,ny)]; zeros(ny,n)]
  b_tilda = [[a12 b1 zeros(nx,ny)]; [zeros(ny,ny+np) eye(ny)]]
  c_tilda = c

  q_tilda = [[q11 zeros(nx,ny)]; zeros(ny,n)]
  u_tilda = [[q12 u1 a21']; [((-1/beta)*a0) zeros(ny,ny+np)]]
  r_tilda = [[q22 u2 a22']; [u2' r b2']; [a22 b2 zeros(ny,ny)]]

  v = q+(a+a')/n

  len = Inf
  iters = 0
  while len > tol && iters < maxiters

    f  = -(r_tilda+beta*b_tilda'v*b_tilda)\(u_tilda'+beta*b_tilda'v*a_tilda)
    vn = q_tilda+u_tilda*f+f'u_tilda'+f'r_tilda*f+beta*(a_tilda+b_tilda*f)'v*(a_tilda+b_tilda*f)

    len   = maximum(abs,vn-v)
    v     = vn

    iters += 1

  end

  if iters <= maxiters
    retcode = true
  else
    retcode = false
  end

  q_orig = q_tilda
  u_orig = [[q12 u1]; zeros(ny,ny+np)]
  r_orig = [[q22 u2]; [u2' r]]

  p = a_tilda+b_tilda*f
  v = dlyap(sqrt(beta)*p',q_orig+u_orig[:,1:(ny+np)]*f[1:(ny+np),:]+f[1:(ny+np),:]'u_orig[:,1:(ny+np)]'+f[1:(ny+np),:]'r_orig[1:(ny+np),1:(ny+np)]*f[1:(ny+np),:])
  k = c_tilda
  h = [[eye(nx) zeros(nx,ny)]; f[1:ny,:]]
  f = f[(ny+1):(ny+np),:]

  # Now re-express the solution by eliminating the multipliers

  p11 = p[1:nx,1:nx]
  p12 = p[1:nx,(nx+1):n]
  p21 = p[(nx+1):n,1:nx]
  p22 = p[(nx+1):n,(nx+1):n]

  k11 = k[1:nx,:]

  h21 = h[(nx+1):n,1:nx]
  h22 = h[(nx+1):n,(nx+1):n]

  f1 = f[1:np,1:nx]
  f2 = f[1:np,(nx+1):n]

  l = [h22; f2]
  g = [h21; f1]

  invl = pinv(l)

  p0                             = eye(2*nx+ny+np)
  p0[1:nx,(2*nx+1):(2*nx+ny+np)] = -p12*invl

  p1                                              = zeros(2*nx+ny+np,2*nx+ny+np)
  p1[1:nx,1:nx]                                   = p11-p12*invl*g
  p1[(nx+1):(2*nx),1:nx]                          = eye(nx)
  p1[(2*nx+1):(2*nx+ny+np),1:nx]                  = g
  p1[(2*nx+1):(2*nx+ny+np),(nx+1):(2*nx)]         = l*(p21-p22*invl*g)
  p1[(2*nx+1):(2*nx+ny+np),(2*nx+1):(2*nx+ny+np)] = l*p22*invl

  p = p0\p1

  h = p[(nx+1):(2*nx+ny),:]
  f = p[(2*nx+ny+1):(2*nx+ny+np),:]
  k = [k11; zeros((nx+ny+np),ns)]

  soln = State_Space_Soln(p,k,h,f,v,sigma,retcode)

  return soln

end
