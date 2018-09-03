function solve_quasi(model::State_Space_Form{T}, obj::State_Space_Objective{T}, commit_prob::T, tol::T, maxiters::S) where {T <: AbstractFloat, S <: Int}

  nx    = copy(model.nx)
  ny    = copy(model.ny)
  beta  = copy(obj.beta)
  sigma = copy(model.sigma)

  a  = copy(model.a)
  b  = copy(model.b)
  c  = copy(model.c)

  q = copy(obj.q)
  u = copy(obj.u)
  r = copy(obj.r)

  n  = nx+ny
  np = size(b,2)

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
  c_tilda = [c; zeros(ny,size(c,2))]

  v       = q+(a+a')/n
  h_reopt = zeros(ny,nx)

  s = [eye(nx) zeros(nx,ny)]  # x_{t} = s*z_{t}

  q_tilda = [[q11 zeros(nx,ny)]; zeros(ny,n)]

  len = Inf
  iters = 0
  while len > tol && iters < maxiters

    u_tilda = [[q12 u1 ((a21-(1-commit_prob)*h_reopt*a11))']; [-eye(ny)/beta zeros(ny,ny+np)]]
    r_tilda = [[q22 u2 ((a22-(1-commit_prob)*h_reopt*a12))']; [u2' r ((b2-(1-commit_prob)*h_reopt*b1))']; [(a22-(1-commit_prob)*h_reopt*a12) (b2-(1-commit_prob)*h_reopt*b1) zeros(ny,ny)]]

    f  = -(r_tilda+beta*b_tilda'*(commit_prob*v+(1-commit_prob)*(s'*pinv(s')*v*pinv(s)*s))*b_tilda)\(u_tilda'+beta*b_tilda'*(commit_prob*v+(1-commit_prob)*(s'*pinv(s')*v*pinv(s)*s))*a_tilda)
    vn = q_tilda+u_tilda*f+f'*u_tilda'+f'*r_tilda*f+beta*(a_tilda+b_tilda*f)'*(commit_prob*v+(1-commit_prob)*(s'*pinv(s')*v*pinv(s)*s))*(a_tilda+b_tilda*f)

    len = maximum(abs,vn-v)
    v = vn

    h_reopt = f[1:ny,1:nx]

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
  k = c
  h = [[eye(nx) zeros(nx,ny)]; f[1:ny,:]]
  f = f[(ny+1):(ny+np),:]

  (v,retcode_doubling) = doubling(sqrt(beta)*p',sqrt(beta)*p,h'q*h+h'u*f+f'u'h+f'r*f,tol,maxiters)

  # v = dlyap(sqrt(beta)*p',h'q*h+h'u*f+f'u'h+f'r*f)

  soln = State_Space_Soln(p,k,h,f,v,sigma,retcode)

  return soln

end

function solve_quasi(model::Generalized_State_Space_Form{T}, obj::State_Space_Objective{T}, commit_prob::T, tol::T, maxiters::S) where {T <: AbstractFloat, S <: Int}

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
  c_tilda = [c; zeros(ny,size(c,2))]

  v       = q+(a+a')/n
  h_reopt = zeros(ny,nx)

  s = [eye(nx) zeros(nx,ny)]  # x_{t} = s*z_{t}

  q_tilda = [[q11 zeros(nx,ny)]; zeros(ny,n)]

  len = Inf
  iters = 0
  while len > tol && iters < maxiters

    u_tilda = [[q12 u1 ((a21-(1-commit_prob)*a0*h_reopt*a11))']; [-a0/beta zeros(ny,ny+np)]]
    r_tilda = [[q22 u2 ((a22-(1-commit_prob)*a0*h_reopt*a12))']; [u2' r ((b2-(1-commit_prob)*a0*h_reopt*b1))']; [(a22-(1-commit_prob)*a0*h_reopt*a12) (b2-(1-commit_prob)*a0*h_reopt*b1) zeros(ny,ny)]]

    f  = -(r_tilda+beta*b_tilda'*(commit_prob*v+(1-commit_prob)*(s'*pinv(s')*v*pinv(s)*s))*b_tilda)\(u_tilda'+beta*b_tilda'*(commit_prob*v+(1-commit_prob)*(s'*pinv(s')*v*pinv(s)*s))*a_tilda)
    vn = q_tilda+u_tilda*f+f'*u_tilda'+f'*r_tilda*f+beta*(a_tilda+b_tilda*f)'*(commit_prob*v+(1-commit_prob)*(s'*pinv(s')*v*pinv(s)*s))*(a_tilda+b_tilda*f)

    len = maximum(abs,vn-v)
    v = vn

    h_reopt = f[1:ny,1:nx]

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
  k = c
  h = [[eye(nx) zeros(nx,ny)]; f[1:ny,:]]
  f = f[(ny+1):(ny+np),:]

  (v,retcode_doubling) = doubling(sqrt(beta)*p',sqrt(beta)*p,h'q*h+h'u*f+f'u'h+f'r*f,tol,maxiters)

  # v = dlyap(sqrt(beta)*p',h'q*h+h'u*f+f'u'h+f'r*f)

  soln = State_Space_Soln(p,k,h,f,v,sigma,retcode)

  return soln

end
