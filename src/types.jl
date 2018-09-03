# Types used for rational expectations

mutable struct Blanchard_Kahn_Form{T<:AbstractFloat,S<:Integer}
  # E_tY[t+1] = A*Y[t] + C*V[t+1]
  nx::S                               # Number of predetermined variables
  ny::S                               # Number of nonpredetermined variables
  a::Array{T,2}                       # Companion matrix
  c::Union{Array{T,2},Array{T,1}}     # Innovation loading matrix
  sigma::Union{Array{T,2},Array{T,1}} # Innovation variance-covariance matrix
end

mutable struct Blanchard_Kahn_Soln{T<:AbstractFloat,S<:Integer}
  p::Union{Array{T,2},Array{T,1}}     # Transition matrix for predetermined variables
  k::Union{Array{T,2},Array{T,1}}     # Innovation loading matrix
  h::Union{Array{T,2},Array{T,1}}     # Decision rule matrix linking nonpredetermined variables to predetermined variables
  sigma::Union{Array{T,2},Array{T,1}} # Innovation variance-covariance matrix
  grc::S                              # Number of eigenvalues greater than cutoff
  soln_type::AbstractString                   # "determinate", "indeterminate", or "explosive"
end

mutable struct Klein_Form{T<:AbstractFloat,S<:Integer}
  # B*E_tY[t+1] = A*Y[t] + C*V[t+1]
  nx::S                               # Number of predetermined variables
  ny::S                               # Number of nonpredetermined variables
  a::Array{T,2}                       # Companion matrix; this is B in Klein's paper
  b::Array{T,2}                       # Lead matrix; this is A in Klein's paper
  c::Union{Array{T,2},Array{T,1}}     # Innovation loading matrix
  sigma::Union{Array{T,2},Array{T,1}} # Innovation variance-covariance matrix
end

mutable struct Klein_Soln{T<:AbstractFloat,S<:Integer}
  p::Union{Array{T,2},Array{T,1}}     # Transition matrix for predetermined variables
  k::Union{Array{T,2},Array{T,1}}     # Innovation loading matrix
  h::Union{Array{T,2},Array{T,1}}     # Decision rule matrix linking nonpredetermined variables to predetermined variables
  sigma::Union{Array{T,2},Array{T,1}} # Innovation variance-covariance matrix
  grc::S                              # Number of eigenvalues greater than cutoff
  soln_type::AbstractString                   # "determinate", "indeterminate", or "explosive"
end

mutable struct Binder_Pesaran_Form{T<:AbstractFloat}
  # A*Y[t] = A1*Y[t-1]+ B*E_tY[t+1] + C*V[t]
  a::Array{T,2}                       # Contemporaneous matrix
  a1::Array{T,2}                      # Lag matrix matrix
  b::Array{T,2}                       # Lead matrix
  c::Union{Array{T,2},Array{T,1}}     # Innovation matrix
  sigma::Union{Array{T,2},Array{T,1}} # Innovation variance-covariance matrix
end

mutable struct Binder_Pesaran_Soln{T<:AbstractFloat,S<:Integer}
  p::Union{Array{T,2},Array{T,1}}     # Transition matrix for predetermined and nonpredetermined variables
  k::Union{Array{T,2},Array{T,1}}     # Innovation loading matrix
  sigma::Union{Array{T,2},Array{T,1}} # Innovation variance-covariance matrix
  grc::S                              # Number of eigenvalues greater than cutoff
  soln_type::AbstractString                   # "determinate", "indeterminate", or "explosive"
end

mutable struct Sims_Form{T<:AbstractFloat}
  # gamma0*Y[t] = gamma1*Y[t-1] + C + psi*V[t] + pi*Eta[t]
  gamma0::Array{T,2}                  # Lead matrix
  gamma1::Array{T,2}                  # Companion matrix
  c::Array{T,1}                       # Constant vector
  psi::Array{T,2}                     # Innovation loading matrix
  pi::Union{Array{T,2},Array{T,1}}    # Expectational-errors loading matrix
  sigma::Union{Array{T,2},Array{T,1}} # Innovation variance-covariance matrix
end

mutable struct Sims_Soln{T<:AbstractFloat,S<:Integer}
  # y(t) = g1*y(t-1) + c + impact*z(t) + ywt*inv(I-fmat*inv(L))*fwt*z(t+1)
  g1::Union{Array{T,2},Array{T,1}}     # Transition matrix
  c::Array{T,1}                        # Constant vector
  impact::Union{Array{T,2},Array{T,1}} # Constant vector
  ywt::Union{Array{T,2},Array{T,1}}    # Innovation loading matrix
  fmat::Union{Array{T,2},Array{T,1}}   # Expectational-errors loading matrix
  fwt::Union{Array{T,2},Array{T,1}}    # Expectational-errors loading matrix
  sigma::Union{Array{T,2},Array{T,1}}  # Innovation variance-covariance matrix
  grc::S                               # Number of eigenvalues greater than cutoff
  iid_exist::Bool                      # true or false
  general_exist::Bool                  # true or false
  uniqueness::Bool                     # true or false
end

mutable struct Gomme_Klein_Form{T<:AbstractFloat,S<:Integer}
  nx::S                                # Number of predetermined variables
  ny::S                                # Number of nonpredetermined variables
  derivs::Array{T,2}                   # Matrix of first derivatives (n*2n, where n = nx+ny)
  hessians::Array{T,2}                 # Matrix of stacked second derivatives ((2n)^2*2n, where n = nx+ny)
  eta::Union{Array{T,2},Array{T,1}}    # Innovation loading-matrix in predetermined equation-block
  sigma::Union{Array{T,2},Array{T,1}} # Innovation variance-covariance matrix
end

mutable struct Gomme_Klein_Soln{T<:AbstractFloat,S<:Integer}
  # x(t+1) = ssh + hx*x(t) + [kron(I,x(t))]'hxx*[kron(I,x(t))] + v(t+1)
  #   y(t) = ssg + gx*x(t) + [kron(I,x(t))]'gxx*[kron(I,x(t))]
  ssh::Array{T,1}                     # Intercepts in predetermined block
  hx::Union{Array{T,2},Array{T,1}}    # Linear component in predetermined block
  hxx::Array{T,2}                     # Quadratic component in predetermined block
  ssg::Array{T,1}                     # Intercepts in predetermined block
  gx::Union{Array{T,2},Array{T,1}}    # Linear component in predetermined block
  gxx::Array{T,2}                     # Quadratic component in predetermined block
  eta::Union{Array{T,2},Array{T,1}}   # Innovation loading-matrix in predetermined equation-block
  sigma::Union{Array{T,2},Array{T,1}} # Innovation variance-covariance matrix
  grc::S                              # Number of eigenvalues greater than cutoff
  soln_type::AbstractString                   # "determinate", "indeterminate", or "explosive"
end

mutable struct Lombardo_Sutherland_Form{T<:AbstractFloat,S<:Integer}
  nx::S                                # Number of predetermined variables
  ny::S                                # Number of nonpredetermined variables
  derivs::Array{T,2}                   # Matrix of first derivatives (n*2n, where n = nx+ny)
  hessians::Array{T,2}                 # Matrix of stacked second derivatives ((2n)^2*2n, where n = nx+ny)
  eta::Union{Array{T,2},Array{T,1}}    # Innovation loading-matrix in predetermined equation-block
  sigma::Union{Array{T,2},Array{T,1}}  # Innovation variance-covariance matrix
end

mutable struct Lombardo_Sutherland_Soln{T<:AbstractFloat,S<:Integer}
  ssh::Union{Array{T,2},Array{T,1}}   # Intercepts in predetermined block
  hx::Union{Array{T,2},Array{T,1}}    # Linear component in predetermined block
  hxx::Array{T,2}                     # Quadratic component in predetermined block
  ssg::Union{Array{T,2},Array{T,1}}   # Intercepts in predetermined block
  gx::Union{Array{T,2},Array{T,1}}    # Linear component in predetermined block
  gxx::Array{T,2}                     # Quadratic component in predetermined block
  eta::Union{Array{T,2},Array{T,1}}   # Innovation loading-matrix in predetermined equation-block
  phi::Array{T,2}
  gamma::Array{T,2}
  psi::Array{T,2}
  sigma::Union{Array{T,2},Array{T,1}} # Innovation variance-covariance matrix
  grc::S                              # Number of eigenvalues greater than cutoff
  soln_type::AbstractString                   # "determinate", "indeterminate", or "explosive"
end

# Types used for optimal policy

mutable struct State_Space_Objective{T<:AbstractFloat}
  beta::T
  q::Array{T,2}
  u::Array{T,2}
  r::Array{T,2}
end

mutable struct State_Space_Form{T<:AbstractFloat,S<:Integer}
  # E_tY[t+1] = A*Y[t] + B*U[t] + C*V[t+1]
  nx::S             # Number of predetermined variables
  ny::S             # Number of nonpredetermined variables
  a::Array{T,2}     # Companion matrix
  b::Array{T,2}     # Policy loadings
  c::Array{T,2}     # Innovation loading matrix
  sigma::Array{T,2} # Innovation variance-covariance matrix
end

mutable struct Generalized_State_Space_Form{T<:AbstractFloat,S<:Integer}
  # E_tY[t+1] = A*Y[t] + B*U[t] + C*V[t+1]
  nx::S             # Number of predetermined variables
  ny::S             # Number of nonpredetermined variables
  a0::Array{T,2}    # Companion matrix
  a::Array{T,2}     # Companion matrix
  b::Array{T,2}     # Policy loadings
  c::Array{T,2}     # Innovation loading matrix
  sigma::Array{T,2} # Innovation variance-covariance matrix
end

mutable struct State_Space_Soln{T<:AbstractFloat}
  p::Array{T,2}     # Transition matrix for predetermined variables
  k::Array{T,2}     # Innovation loading matrix
  h::Array{T,2}
  f::Array{T,2}     # Decision rule matrix linking nonpredetermined variables to predetermined variables
  v::Array{T,2}
  sigma::Array{T,2} # Innovation variance-covariance matrix
  converged::Bool
end

mutable struct Structural_Objective{T<:AbstractFloat}
  beta::T
  q::Array{T,2}
  r::Array{T,2}
end

mutable struct Structural_Form{T<:AbstractFloat}
  # E_tY[t+1] = A*Y[t] + B*U[t] + C*V[t+1]
  a0::Array{T,2}     # Companion matrix
  a1::Array{T,2}     # Policy loadings
  a2::Array{T,2}     # Innovation loading matrix
  a3::Array{T,2}     # Policy loadings
  a5::Array{T,2}     # Innovation loading matrix
  sigma::Array{T,2}  # Innovation variance-covariance matrix
end

mutable struct Generalized_Structural_Form{T<:AbstractFloat}
  # E_tY[t+1] = A*Y[t] + B*U[t] + C*V[t+1]
  a0::Array{T,2}     # Companion matrix
  a1::Array{T,2}     # Policy loadings
  a2::Array{T,2}     # Innovation loading matrix
  a3::Array{T,2}     # Policy loadings
  a4::Array{T,2}     # Policy loadings
  a5::Array{T,2}     # Innovation loading matrix
  sigma::Array{T,2}  # Innovation variance-covariance matrix
end

mutable struct Structural_Soln{T<:AbstractFloat}
  p::Array{T,2}     # Transition matrix for predetermined variables
  k::Array{T,2}     # Innovation loading matrix
  v::Array{T,2}
  sigma::Array{T,2} # Innovation variance-covariance matrix
  converged::Bool
end

# Types used for impulse response functions

Perturbable_Soln = Union{Blanchard_Kahn_Soln,Klein_Soln,Binder_Pesaran_Soln,Gomme_Klein_Soln,Lombardo_Sutherland_Soln,State_Space_Soln,Structural_Soln}
Second_Order_State_Space_Soln = Union{Gomme_Klein_Soln,Lombardo_Sutherland_Soln}
