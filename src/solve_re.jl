function solve_re(model::Blanchard_Kahn_Form{T}, cutoff::T) where T <: AbstractFloat

  nx    = copy(model.nx)
  ny    = copy(model.ny)
  a     = copy(model.a)
  c     = copy(model.c)
  sigma = copy(model.sigma)

  n = nx+ny

  r = schur(complex(a))

  # Reorder the eigenvalues so that those with modulus greater than "cutoff" reside at the bottom.

  sel = (abs.(r.values) .< cutoff)
  ordschur!(r,sel)
  s = r.T
  q = r.Z' # So now q*a*q' = s

  # Calculate the number of unstable eigenvalues

  grc = n-sum(sel)

  soln_type = "determinate"
  if grc < ny
    soln_type = "indeterminate"
  elseif grc > ny
    soln_type = "explosive"
  end

  # Construct the rational expectations equilibrium by eliminating the unstable dynamics

  q11 = q[1:nx,1:nx]
  q12 = q[1:nx,(nx+1):n]
  q21 = q[(nx+1):n,1:nx]
  q22 = q[(nx+1):n,(nx+1):n]
  s11 = s[1:nx,1:nx]

  p = real((q11-(q12/q22)*q21)\s11*(q11-(q12/q22)*q21))
  k = c[1:nx,:]
  h = -real(q22\q21)

  soln = Blanchard_Kahn_Soln(p,k,h,sigma,grc,soln_type)

  return soln

end

function solve_re(model::Blanchard_Kahn_Form{T}, cutoff::T, tol::T) where T <: AbstractFloat

  nx    = copy(model.nx)
  ny    = copy(model.ny)
  a     = copy(model.a)
  c     = copy(model.c)
  sigma = copy(model.sigma)

  n = nx+ny

  s = eigvals(a)

  grc = 0
  for i = 1:n

    if norm(s[i]) > cutoff
      grc += 1
    end

  end

  soln_type = "determinate"
  if grc < ny
    soln_type = "indeterminate"
  elseif grc > ny
    soln_type = "explosive"
  end

  # Construct the rational expectations equilibrium

  a11 = a[1:nx,1:nx]
  a12 = a[1:nx,(nx+1):n]
  a21 = a[(nx+1):n,1:nx]
  a22 = a[(nx+1):n,(nx+1):n]

  h0 = -a22\a21
  h = Array{T}(undef,ny,nx)
  len = Inf
  while len > tol

    h = (a22-h0*a12)\(h0*a11-a21)
    len = maximum(abs,h-h0)
    h0 = h

  end

  p = a11+a12*h
  k = c[1:nx,:]

  soln = Blanchard_Kahn_Soln(p,k,h,sigma,grc,soln_type)

  return soln

end

function solve_re(model::Klein_Form{T}, cutoff::T) where T <: AbstractFloat

  nx    = copy(model.nx)
  ny    = copy(model.ny)
  a     = copy(model.a)
  b     = copy(model.b)
  c     = copy(model.c)
  sigma = copy(model.sigma)

  n = nx+ny

  r = schur(complex(a),complex(b))

  # Reorder the generalized eigenvalues so that those with modulus greater than "cutoff" reside at the bottom.

  sel = (abs.(diag(r.S)./diag(r.T)).<cutoff)
  ordschur!(r,sel)
  s = r.S
  t = r.T
  q = r.Q'  # So now q*a*z = s and q*b*z = t
  z = r.Z

  # Calculate the number of unstable eigenvalues

  grc = n-sum(sel)

  soln_type = "determinate"
  if grc < ny
    soln_type = "indeterminate"
  elseif grc > ny
    soln_type = "explosive"
  end

  # Construct the rational expectations equilibrium by eliminating the unstable dynamics

  z11 = z[1:nx,1:nx]
  z21 = z[(nx+1):n,1:nx]
  t11 = t[1:nx,1:nx]
  s11 = s[1:nx,1:nx]

  p = real((z11/t11)*(s11/z11))
  k = c[1:nx,:]
  h = real(z21/z11)

  soln = Klein_Soln(p,k,h,sigma,grc,soln_type)

  return soln

end

function solve_re(model::Binder_Pesaran_Form{T}, cutoff::T) where T <: AbstractFloat

  a     = copy(model.a)
  a1    = copy(model.a1)
  b     = copy(model.b)
  c     = copy(model.c)
  sigma = copy(model.sigma)

  nx = size(a1,1)
  ny = size(b,1)

  n = nx+ny

  a = [zeros(nx,nx) I; -a1 a]
  b = [I zeros(nx,ny); zeros(ny,nx) b]
  c = [zeros(nx,size(c,2)); -c]

  r = schur(complex(a),complex(b))

  # Reorder the generalized eigenvalues so that those with modulus greater than "cutoff" reside at the bottom.

  sel = (abs.(diag(r.S)./diag(r.T)).<cutoff)
  ordschur!(r,sel)
  s = r.S
  t = r.T
  q = r.Q'  # So now q*a*z = s and q*b*z = t
  z = r.Z

  # Calculate the number of unstable eigenvalues

  grc = n-sum(sel)

  soln_type = "determinate"
  if grc < ny
    soln_type = "indeterminate"
  elseif grc > ny
    soln_type = "explosive"
  end

  # Construct the rational expectations equilibrium by eliminating the unstable dynamics

  z11 = z[1:nx,1:nx]
  z12 = z[1:nx,(nx+1):n]
  z21 = z[(nx+1):n,1:nx]
  z22 = z[(nx+1):n,(nx+1):n]
  t11 = t[1:nx,1:nx]
  s11 = s[1:nx,1:nx]
  s12 = s[1:nx,(nx+1):n]
  s22 = s[(nx+1):n,(nx+1):n]

  qc  = -q*c
  qcs = qc[1:nx,:]       # The part of qc that is in the stable block
  qcu = qc[(nx+1):n,:]   # The part of qc that is in the unstable block

  p = real((z11/t11)*(s11/z11))
  k = real((z11/t11)*(s12/s22)*qcu-(z11/t11)*(s11/z11)*(z12/s22)*qcu-(z11/t11)*qcs)

  soln = Binder_Pesaran_Soln(p,k,sigma,grc,soln_type)

  return soln

end

function solve_re(model::Binder_Pesaran_Form{T}, cutoff::T, tol::T) where T <: AbstractFloat

  a     = copy(model.a)
  a1    = copy(model.a1)
  b     = copy(model.b)
  c     = copy(model.c)
  sigma = copy(model.sigma)

  p0 = copy(a1)
  p  = copy(p0)

  len = Inf
  while len > tol

    p = (a-b*p0)\a1
    len = maximum(abs,p-p0)
    p0 = p

  end

  k = (a-b*p0)\c

  nx = size(a1,1)
  ny = size(b,1)
  aa = [zeros(nx,nx) I; -a1 a]
  bb = [I zeros(nx,ny); zeros(ny,nx) b]

  lambda = eigvals(bb,aa)

  grc = 0
  for i = 1:(nx+ny)

    if norm(lambda[i]) > cutoff
      grc += 1
    end

  end

  soln_type = "determinate"
  if grc < ny
    soln_type = "indeterminate"
  elseif grc > ny
    soln_type = "explosive"
  end

  soln = Binder_Pesaran_Soln(p,k,sigma,grc,soln_type)

  return soln

end

function solve_re(model::Sims_Form{T}, cutoff::T) where T <: AbstractFloat

  # This function is an implementation of a routine written by
  # Christopher A. Sims (http://www.princeton.edu/~sims/)

  # Permission to release this function under a MIT licence has
  # been granted by Chris Sims, 23/8/2015.

  gamma0 = copy(model.gamma0)
  gamma1 = copy(model.gamma1)
  c      = copy(model.c)
  psi    = copy(model.psi)
  pi     = copy(model.pi)
  sigma  = copy(model.sigma)

  n  = size(gamma0,1)

  r = schur(complex(gamma1),complex(gamma0))

  # Reorder the generalized eigenvalues so that those with modulus greater than "cutoff" reside at the bottom.

  sel = (abs.(diag(r.S)./diag(r.T)).<cutoff)
  ordschur!(r,sel)
  s = r.S
  t = r.T
  q = r.Q'  # So now q*a*z = s and q*b*z = t
  z = r.Z

  # Calculate the number of unstable eigenvalues

  grc = n-sum(sel)

  # Define some ingredients that are needed for any solution

  q1 = q[1:(n-grc),:]
  q2 = q[(n-grc+1):n,:]
  z1 = z[:,1:(n-grc)]'
  z2 = z[:,(n-grc+1):n]'
  t2 = t[(n-grc+1):n,(n-grc+1):n]
  s2 = s[(n-grc+1):n,(n-grc+1):n]
  q2pi  = q2*pi
  q2psi = q2*psi

  # Check for the existance of a solution

  (pi_rows,pi_cols)   = size(pi)
  (psi_rows,psi_cols) = size(psi)

  if grc == 0
    q2pi    = zeros(0,pi_cols)
    uq2pi   = zeros(0,0)
    dq2pi   = zeros(0,0)
    vq2pi   = zeros(pi_cols,0)

    q2psi   = zeros(0,psi_cols)
    uq2psi  = zeros(0,0)
    dq2psi  = zeros(0,0)
    vq2psi  = zeros(psi_cols,0)

    non_zero_index_psi = []
  else
    (uq2pi,dq2pi,vq2pi) = svd(q2pi)
    non_zero_index_pi = findall(dq2pi .> 2*eps(T)*norm(q2pi))
    uq2pi = uq2pi[:,non_zero_index_pi]
    vq2pi = vq2pi[:,non_zero_index_pi]
    dq2pi = Matrix(Diagonal(dq2pi[non_zero_index_pi]))

    (uq2psi,dq2psi,vq2psi) = svd(q2psi)
    non_zero_index_psi = findall(dq2psi .> 2*eps(T)*norm(q2psi))
    uq2psi = uq2psi[:,non_zero_index_psi]
    vq2psi = vq2psi[:,non_zero_index_psi]
    dq2psi = Matrix(Diagonal(dq2psi[non_zero_index_psi]))
  end

  if isempty(non_zero_index_psi)
    iid_exist = true
  else
    iid_exist = (norm(uq2psi - uq2pi*uq2pi'*uq2psi) < 2*eps(T)*norm(uq2psi))
  end

  if iid_exist == false
    general_exist = false
  else
    t2iq2psi = t2\q2psi
    p        = t2\s2
    h        = t2iq2psi
    for i = 1:n

      h = [p*h t2iq2psi]

    end
    (uh,dh,vh) = svd(h)
    non_zero_index_dh = findall(dh .> 2*eps(T)*norm(h))
    uh = uh[:,non_zero_index_dh]
    vh = vh[:,non_zero_index_dh]
    dh = diagm(dh[non_zero_index_dh])
    general_exist = (norm(uh-uq2pi*uq2pi'*uh) < 2*eps(T)*norm(uh,1))
  end

  # Check for uniqueness

  if grc == n
    q1pi  = zeros(0,pi_cols)
    uq1pi = zeros(0,0)
    vq1pi = zeros(pi_cols,0)
    dq1pi = zeros(0,0)
  else
    q1pi = q1*pi
    (uq1pi,dq1pi,vq1pi) = svd(q1pi)
    non_zero_index_dq1pi = findall(dq1pi .> 2*eps(T)*norm(q1pi))
    uq1pi = uq1pi[:,non_zero_index_dq1pi]
    vq1pi = vq1pi[:,non_zero_index_dq1pi]
    dq1pi = Matrix(Diagonal(dq1pi[non_zero_index_dq1pi]))
  end

  if isempty(vq1pi)
    unique = true
  else
    resid = vq1pi-vq2pi*vq2pi'vq1pi
    (u1,d1,v1) = svd(resid)
    unique = ((sum(abs,d1) > 2*eps(T)*n) != true)
  end
  tmat = [I -(uq2pi*(dq2pi\vq2pi')*vq1pi*dq1pi*uq1pi')']
  g0   = [tmat*t; zeros(grc,(n-grc)) I]
  g1   = [tmat*s; zeros(grc,n)]

  g0i    = inv(g0)
  g1     = real(z*g0i*g1*z')
  c      = real(z*g0i*[tmat*q*c; (t[(n-grc+1):n,(n-grc+1):n]-s[(n-grc+1):n,(n-grc+1):n])\q2*c])
  impact = real(z*g0i*[tmat*q*psi; zeros(grc,size(psi,2))])
  fmat   = real(s[(n-grc+1):n,(n-grc+1):n]\t[(n-grc+1):n,(n-grc+1):n])
  fwt    = -real(s[(n-grc+1):n,(n-grc+1):n]\q2*psi)
  ywt    = real(z*g0i[:,(n-grc+1):n])

  soln = Sims_Soln(g1,c,impact,ywt,fmat,fwt,sigma,grc,iid_exist,general_exist,unique)

  return soln

end

function solve_re(model::Gomme_Klein_Form, cutoff::T) where T <: AbstractFloat

  nx     = copy(model.nx)
  ny     = copy(model.ny)
  deriv1 = copy(model.derivs)
  deriv2 = copy(model.hessians)
  eta    = copy(model.eta)
  sigma  = copy(model.sigma)

  n = nx+ny

  # Reorder the second derivatives to put them in the order that Gomme and Klein use
  # Probably do away with this in a subsequent update

  hes = [deriv2[(n+1):2*n,:]; deriv2[1:n,:]]
  for i = 2:n
    hes = [hes; deriv2[2*n*(i-1)+n+1:2*n*i,:]; deriv2[2*n*(i-1)+1:2*n*(i-1)+n,:]]
  end
  hes = [hes[:,(n+1):2*n] hes[:,1:n]]

  # Construct partitioned first derivative matrices

  fx  = deriv1[1:n,1:nx]
  fy  = deriv1[1:n,(nx+1):n]
  fxp = deriv1[1:n,(n+1):(n+nx)]
  fyp = deriv1[1:n,(n+nx+1):(2*n)]

  # Solve the linearized first-order system

  a =  [fx fy]
  b = -[fxp fyp]
  c =  [eta; zeros(ny,size(eta,2))]

  first_order_model = Klein_Form(nx,ny,a,b,c,model.sigma)
  first_order_soln  = solve_re(first_order_model,cutoff)

  hx = first_order_soln.p
  gx = first_order_soln.h
  grc = first_order_soln.grc

  soln_type = "determinate"
  if grc < ny
    soln_type = "indeterminate"
  elseif grc > ny
    soln_type = "explosive"
  end

  # Set up the LP problem needed to construct the matrices on the second-order states

  m = [hx; gx*hx; I; gx]
  q = kron(Matrix(1.0I,n,n),m')*hes*m
  b1 = kron(fxp,Matrix(1.0I,nx,nx))
  b2 = kron(fyp,Matrix(1.0I,nx,nx))
  b4 = kron(fy,Matrix(1.0I,nx,nx))
  c1 = kron(Matrix(1.0I,ny,ny),hx')
  c2 = kron(gx,Matrix(1.0I,nx,nx))

  qq = [kron(Matrix(1.0I,nx,nx),b4)+kron(hx',b2*c1) kron(Matrix(1.0I,nx,nx),b1+b2*c2)]
  xx = -qq\vec(q)
  gxx = reshape(xx[1:nx^2*ny],nx*ny,nx)
  hxx = reshape(xx[nx^2*ny+1:n*nx*nx],nx^2,nx)

  # Set up the LP problem needed to construct the intercepts, which contain the volatility effects

  qq = [fxp+fyp*gx fyp+fy]
  nn = [I; gx; zeros(n,nx)]
  q = fyp*tracem(kron(Matrix(1.0I,ny,ny),eta*(sigma')*eta')*gxx)+tracem(kron(Matrix(1.0I,n,n),nn')*hes*nn*eta*(sigma')*eta')
  ss = -qq\q
  ssh = ss[1:nx]
  ssg = ss[nx+1:n]

  soln = Gomme_Klein_Soln(ssh,hx,hxx,ssg,gx,gxx,eta,sigma,grc,soln_type)

  return soln

end

function solve_re(model::Lombardo_Sutherland_Form, cutoff::T) where T <: AbstractFloat

  nx     = copy(model.nx)
  ny     = copy(model.ny)
  deriv1 = copy(model.derivs)
  deriv2 = copy(model.hessians)
  eta    = copy(model.eta)
  sigma  = copy(model.sigma)

  n  = nx+ny
  ns = size(eta,2)

  # Solve the first order system

  a1 = -deriv1[1:n,(n+1):2*n]
  a2 = deriv1[1:n,1:n]
  c  = [eta; zeros(ny,size(eta,2))]

  r = schur(complex(a2),complex(a1))

  # Reorder the generalized eigenvalues so that those with modulus greater than "cutoff" reside at the bottom.

  sel = (abs.(diag(r.S)./diag(r.T)).<cutoff)
  ordschur!(r,sel)
  s = r.S
  t = r.T
  q = r.Q'  # So now q*a*z = s and q*b*z = t
  z = r.Z

  # Calculate the number of unstable eigenvalues

  grc = n-sum(sel)

  soln_type = "determinate"
  if grc < ny
    soln_type = "indeterminate"
  elseif grc > ny
    soln_type = "explosive"
  end

  # Construct the rational expectations equilibrium by eliminating the unstable dynamics

  z11 = z[1:nx,1:nx]
  z12 = z[1:nx,nx+1:n]
  z21 = z[(nx+1):n,1:nx]
  z22 = z[(nx+1):n,(nx+1):n]

  t11 = t[1:nx,1:nx]
  t12 = t[1:nx,nx+1:n]
  t22 = t[nx+1:n,nx+1:n]

  s11 = s[1:nx,1:nx]
  s12 = s[1:nx,nx+1:n]
  s22 = s[nx+1:n,nx+1:n]

  phi   = real((z11/t11)*(s11/z11))
  gamma = c[1:nx,:]
  omega = [Matrix(1.0I,nx,nx); real(z21/z11)]  # Defines [x; y] as a function of [x]

  soln_type = "determinate"
  if grc < ny
    soln_type = "indeterminate"
  elseif grc > ny
    soln_type = "explosive"
  end

  # Set up the matrices for the second order terms

  for i = 1:n
    for  j = 1:2*n
      deriv2[(i-1)*2*n+j,j] = deriv2[(i-1)*2*n+j,j]/2
    end
    deriv2[(i-1)*2*n+1:i*2*n,:] = 2*deriv2[(i-1)*2*n+1:i*2*n,:]
  end

  v2vh_n = vec_to_vech(n)

  a4 = v2vh_n*vec(deriv2[1:n,1:n])
  a5 = v2vh_n*vec(deriv2[(n+1):2*n,(n+1):2*n])
  for i = 2:n
    a4 = [a4 v2vh_n*vec(deriv2[(i-1)*2*n+1:(i-1)*2*n+n,1:n])]
    a5 = [a5 v2vh_n*vec(deriv2[(i-1)*2*n+n+1:i*2*n,(n+1):2*n])]
  end

  a4 = a4'
  a5 = a5'

  vh2v_nx = vech_to_vec(nx)
  v2vh_nx = vec_to_vech(nx)
  vh2v_ns = vech_to_vec(ns)
  v2vh_ns = vec_to_vech(ns)
  p       = permutation(nx,ns)

  r           = v2vh_n*kron(omega,omega)*vh2v_nx            # n*(n+1)/2   times nx*(nx+1)/2
  phi_tilda   = v2vh_nx*kron(phi,phi)*vh2v_nx               # nx*(nx+1)/2 times nx*(nx+1)/2
  gamma_tilda = v2vh_nx*kron(gamma,gamma)*vh2v_ns           # nx*(nx+1)/2 times ns*(ns+1)/2
  psi_tilda   = v2vh_nx*(kron(phi,gamma)+kron(gamma,phi)*p) # nx*(nx+1)/2 times nx*ns

  g = a4*r+a5*r*phi_tilda                                   # n times nx*(nx+1)/2
  h = a5*r*gamma_tilda                                      # n times ns*(ns+1)/2

  sig = v2vh_ns*vec(sigma)                                  # ns*(ns+1)/2 times ns*ns

  # The next step is to resolve the system using the q and z matrices obtained earlier

  g = q*g
  h = q*h

  gs = g[1:nx,:]
  gu = g[nx+1:n,:]

  hs = h[1:nx,:]
  hu = h[nx+1:n,:]

  m = reshape((kron(phi_tilda',t22)-kron(Matrix(1.0I,round(Int,(nx*(nx+1)/2)),round(Int,(nx*(nx+1)/2))),s22))\vec(gu),ny,round(Int,(nx*(nx+1)/2)))

  hx  = real(z11*(t11\s11)/z11)
  hxx = real(-(z11*(t11\s11)/z11)*z12*m + z11*(t11\(s12*m-t12*m*phi_tilda+gs)) + z12*m*phi_tilda)

  gx  = real(z21/z11)
  gxx = real((z22-z21*(z11\z12))*m)

  vbar = (Matrix(1.0I,round(Int,(nx*(nx+1)/2)),round(Int,(nx*(nx+1)/2)))-phi_tilda)\gamma_tilda
  m3   = (Matrix(1.0I,ny,ny)-s22\t22)\(s22\(gu*vbar+hu))
  m2   = reshape((Matrix(1.0I,round(Int,(nx*(nx+1)/2)*ny),round(Int,(nx*(nx+1)/2)*ny))-kron(phi_tilda',s22\t22))\(vec(s22\gu)),ny,round(Int,(nx*(nx+1)/2)))

  p4 = -real(z22'\(m3-m2*vbar))
  p3 = -z22'\m2
  p2 = -z22'\z12'

  r4 = (t11*z21'+t12*z22')*p4
  r3 = (t11*z21'+t12*z22')*p3
  r2 = t11*z11'+t12*z12'+(t11*z21'+t12*z22')*p2

  d4 = (s11*z21'+s12*z22')*p4+hs

  ssg = p4*sig
  ssh = real(r2\(d4-r4-r3*gamma_tilda))*sig

  soln = Lombardo_Sutherland_Soln(ssh,hx,hxx,ssg,gx,gxx,eta,phi_tilda,gamma_tilda,psi_tilda,sigma,grc,soln_type)

  return soln

end
