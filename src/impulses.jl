function impulses(solution::Perturbable_Soln, impulse_length::S, innovation_to_shock::S) where S <: Int

  if isa(solution,Blanchard_Kahn_Soln) == true
    responses = compute_first_order_state_space_re_impulses(solution,impulse_length,innovation_to_shock)
  elseif isa(solution,Klein_Soln) == true
      responses = compute_first_order_state_space_re_impulses(solution,impulse_length,innovation_to_shock)
  elseif isa(solution,State_Space_Soln) == true
    responses = compute_state_space_op_impulses(solution,impulse_length,innovation_to_shock)
  elseif isa(solution,Binder_Pesaran_Soln) == true
    responses = compute_structural_re_impulses(solution,impulse_length,innovation_to_shock)
  elseif isa(solution,Structural_Soln) == true
    responses = compute_structural_op_impulses(solution,impulse_length,innovation_to_shock)
  elseif isa(solution,Second_Order_State_Space_Soln) == true
    if isa(solution,Lombardo_Sutherland_Soln) == true
      solution = convert_second_order(solution)
    end
    responses = compute_second_order_state_space_re_impulses(solution,impulse_length,innovation_to_shock)
  end

end

function compute_first_order_state_space_re_impulses(solution,impulse_length,innovation_to_shock)

  p     = copy(solution.p)
  k     = copy(solution.k)
  h     = copy(solution.h)
  sigma = copy(solution.sigma)

  nx = size(p,1)
  ny = size(h,1)

  s = cholesky(sigma[:,:]).U'

  responses = zeros(nx+ny,impulse_length)
  responses[1:nx,1]           = (k*s)[:,innovation_to_shock]
  responses[(nx+1):(nx+ny),1] = h*responses[1:nx,1]

  for i = 2:impulse_length

    responses[1:nx,i]           = p*responses[1:nx,i-1]
    responses[(nx+1):(nx+ny),i] = h*responses[1:nx,i]

  end

  return responses

end

function compute_state_space_op_impulses(solution,impulse_length,innovation_to_shock)

  p     = copy(solution.p)
  k     = copy(solution.k)
  h     = copy(solution.h)
  f     = copy(solution.f)
  sigma = copy(solution.sigma)

  nx = size(p,1)
  ny = size(h,1)
  np = size(f,1)

  s = chol(sigma[:,:])'

  responses = zeros(nx+ny+np,impulse_length)
  responses[1:nx,1]                 = (k*s)[:,innovation_to_shock]
  responses[(nx+1):(nx+ny),1]       = h*responses[1:nx,1]
  responses[(nx+ny+1):(nx+ny+np),1] = f*responses[1:nx,1]

  for i = 2:impulse_length

    responses[1:nx,i]                 = p*responses[1:nx,i-1]
    responses[(nx+1):(nx+ny),i]       = h*responses[1:nx,i]
    responses[(nx+ny+1):(nx+ny+np),i] = f*responses[1:nx,i]

  end

  return responses

end

function compute_structural_re_impulses(solution,impulse_length,innovation_to_shock)

  p     = copy(solution.p)
  k     = copy(solution.k)
  sigma = copy(solution.sigma)

  s = cholesky(sigma[:,:]).U'

  n = size(p,1)

  responses = zeros(n,impulse_length)
  responses[1:n,1] = (k*s)[:,innovation_to_shock]

  for i = 2:impulse_length

    responses[1:n,i] = p*responses[1:n,i-1]

  end

  return responses

end

function compute_structural_op_impulses(solution,impulse_length,innovation_to_shock)

  p     = copy(solution.p)
  k     = copy(solution.k)
  sigma = copy(solution.sigma)

  s = chol(sigma[:,:])'

  n = size(p,1)

  responses = zeros(n,impulse_length)
  responses[1:n,1] = (k*s)[:,innovation_to_shock]

  for i = 2:impulse_length

    responses[1:n,i] = p*responses[1:n,i-1]

  end

  return responses

end

function compute_second_order_state_space_re_impulses(solution,impulse_length,innovation_to_shock)

  hx    = copy(solution.hx)
  gx    = copy(solution.gx)
  hxx   = copy(solution.hxx)
  gxx   = copy(solution.gxx)
  eta   = copy(solution.eta)
  sigma = copy(solution.sigma)

  s = cholesky(sigma[:,:]).U'

  nx = size(hx,1)
  ny = size(gx,1)

  responses = zeros(nx+ny,impulse_length)
  responses[1:nx,1]           = (eta*convert(Array,s))[:,innovation_to_shock]
  responses[(nx+1):(nx+ny),1] = gx*responses[1:nx,1] + 0.5*(kron(Matrix(1.0I,ny,ny),responses[1:nx,1]'))*gxx*responses[1:nx,1]

  for i = 2:impulse_length

    responses[1:nx,i]           = hx*responses[1:nx,i-1] + 0.5*(kron(Matrix(1.0I,nx,nx),responses[1:nx,i-1]'))*hxx*responses[1:nx,i-1]
    responses[(nx+1):(nx+ny),i] = gx*responses[1:nx,i]   + 0.5*(kron(Matrix(1.0I,ny,ny),responses[1:nx,i-1]'))*gxx*responses[1:nx,i-1]

  end

  return responses

end
