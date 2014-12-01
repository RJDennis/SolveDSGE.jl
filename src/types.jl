type Blanchard_Kahn_Form{T<:FloatingPoint,S<:Integer}
  # E_tY[t+1] = A*Y[t] + C*V[t+1]
  nx::S                               # Number of predetermined variables
  ny::S                               # Number of nonpredetermined variables
  a::Array{T,2}                       # Companion matrix
  c::Union(Array{T,2},Array{T,1})     # Innovation loading matrix
  sigma::Union(Array{T,2},Array{T,1}) # Innovation variance-covariance matrix
end

type Blanchard_Kahn_Soln{T<:FloatingPoint,S<:Integer}
  p::Union(Array{T,2},Array{T,1})     # Transition matrix for predetermined variables
  k::Union(Array{T,2},Array{T,1})     # Innovation loading matrix
  f::Union(Array{T,2},Array{T,1})     # Decision rule matrix linking nonpredetermined variables to predetermined variables
  sigma::Union(Array{T,2},Array{T,1}) # Innovation variance-covariance matrix
  grc::S                              # Number of eigenvalues greater than cutoff
  soln_type::String                   # "determinate", "indeterminate", or "explosive"
end

type Klein_Form{T<:FloatingPoint,S<:Integer}
  # B*E_tY[t+1] = A*Y[t] + C*V[t+1]
  nx::S                               # Number of predetermined variables
  ny::S                               # Number of nonpredetermined variables
  a::Array{T,2}                       # Companion matrix; this is B in Klein's paper
  b::Array{T,2}                       # Lead matrix; this is A in Klein's paper
  c::Union(Array{T,2},Array{T,1})     # Innovation loading matrix
  sigma::Union(Array{T,2},Array{T,1}) # Innovation variance-covariance matrix
end

type Klein_Soln{T<:FloatingPoint,S<:Integer}
  p::Union(Array{T,2},Array{T,1})     # Transition matrix for predetermined variables
  k::Union(Array{T,2},Array{T,1})     # Innovation loading matrix
  f::Union(Array{T,2},Array{T,1})     # Decision rule matrix linking nonpredetermined variables to predetermined variables
  sigma::Union(Array{T,2},Array{T,1}) # Innovation variance-covariance matrix
  grc::S                              # Number of eigenvalues greater than cutoff
  soln_type::String                   # "determinate", "indeterminate", or "explosive"
end

type Structural_Form{T<:FloatingPoint}
  # A*Y[t] = A1*Y[t-1]+ B*E_tY[t+1] + C*V[t]
  a::Array{T,2}                       # Contemporaneous matrix
  a1::Array{T,2}                      # Lag matrix matrix
  b::Array{T,2}                       # Lead matrix
  c::Union(Array{T,2},Array{T,1})     # Innovation matrix
  sigma::Union(Array{T,2},Array{T,1}) # Innovation variance-covariance matrix
end

type Structural_Soln{T<:FloatingPoint,S<:Integer}
  p::Union(Array{T,2},Array{T,1})     # Transition matrix for predetermined and nonpredetermined variables
  k::Union(Array{T,2},Array{T,1})     # Innovation loading matrix
  sigma::Union(Array{T,2},Array{T,1}) # Innovation variance-covariance matrix
  grc::S                              # Number of eigenvalues greater than cutoff
  soln_type::String                   # "determinate", "indeterminate", or "explosive"
end

type Sims_Form{T<:FloatingPoint}
  # gamma0*Y[t] = gamma1*Y[t-1] + C + psi*V[t] + pi*Eta[t]
  gamma0::Array{T,2}                  # Lead matrix
  gamma1::Array{T,2}                  # Companion matrix
  c::Array{T,1}                       # Constant vector
  psi::Array{T,2}                     # Innovation loading matrix
  pi::Union(Array{T,2},Array{T,1})    # Expectational-errors loading matrix
  sigma::Union(Array{T,2},Array{T,1}) # Innovation variance-covariance matrix
end

type Sims_Soln{T<:FloatingPoint,S<:Integer}
  # y(t) = g1*y(t-1) + c + impact*z(t) + ywt*inv(I-fmat*inv(L))*fwt*z(t+1)
  g1::Union(Array{T,2},Array{T,1})     # Transition matrix
  c::Array{T,1}                        # Constant vector
  impact::Union(Array{T,2},Array{T,1}) # Constant vector
  ywt::Union(Array{T,2},Array{T,1})    # Innovation loading matrix
  fmat::Union(Array{T,2},Array{T,1})   # Expectational-errors loading matrix
  fwt::Union(Array{T,2},Array{T,1})    # Expectational-errors loading matrix
  sigma::Union(Array{T,2},Array{T,1})  # Innovation variance-covariance matrix
  grc::S                               # Number of eigenvalues greater than cutoff
  iid_exist::Bool                      # true or false
  general_exist::Bool                  # true or false
  uniqueness::Bool                     # true or false
end

type Gomme_Klein_Form{T<:FloatingPoint,S<:Integer}
  nx::S                                # Number of predetermined variables
  ny::S                                # Number of nonpredetermined variables
  derivs::Array{T,2}                   # Matrix of first derivatives (n*2n, where n = nx+ny)
  hessians::Array{T,2}                 # Matrix of stacked second derivatives ((2n)^2*2n, where n = nx+ny)
  eta::Union(Array{T,2},Array{T,1})    # Innovation loading-matrix in predetermined equation-block
  sigma::Union(Array{T,2},Array{T,1}) # Innovation variance-covariance matrix
end

type Gomme_Klein_Soln{T<:FloatingPoint,S<:Integer}
  # x(t+1) = ssh + hx*x(t) + [kron(I,x(t))]'hxx*[kron(I,x(t))] + v(t+1)
  #   y(t) = ssg + gx*x(t) + [kron(I,x(t))]'gxx*[kron(I,x(t))]
  ssh::Array{T,1}                     # Intercepts in predetermined block
  hx::Union(Array{T,2},Array{T,1})    # Linear component in predetermined block
  hxx::Array{T,2}                     # Quadratic component in predetermined block
  ssg::Array{T,1}                     # Intercepts in predetermined block
  gx::Union(Array{T,2},Array{T,1})    # Linear component in predetermined block
  gxx::Array{T,2}                     # Quadratic component in predetermined block
  sigma::Union(Array{T,2},Array{T,1}) # Innovation variance-covariance matrix
  grc::S                              # Number of eigenvalues greater than cutoff
  soln_type::String                   # "determinate", "indeterminate", or "explosive"
end

type Lombardo_Sutherland_Form{T<:FloatingPoint,S<:Integer}
  nx::S                                # Number of predetermined variables
  ny::S                                # Number of nonpredetermined variables
  derivs::Array{T,2}                   # Matrix of first derivatives (n*2n, where n = nx+ny)
  hessians::Array{T,2}                 # Matrix of stacked second derivatives ((2n)^2*2n, where n = nx+ny)
  eta::Union(Array{T,2},Array{T,1})    # Innovation loading-matrix in predetermined equation-block
  sigma::Union(Array{T,2},Array{T,1})  # Innovation variance-covariance matrix
end

type Lombardo_Sutherland_Soln{T<:FloatingPoint,S<:Integer}
  ssh::Union(Array{T,2},Array{T,1})   # Intercepts in predetermined block
  hx::Union(Array{T,2},Array{T,1})    # Linear component in predetermined block
  hxx::Array{T,2}                     # Quadratic component in predetermined block
  ssg::Union(Array{T,2},Array{T,1})   # Intercepts in predetermined block
  gx::Union(Array{T,2},Array{T,1})    # Linear component in predetermined block
  gxx::Array{T,2}                     # Quadratic component in predetermined block
  phi::Array{T,2}
  gamma::Array{T,2}
  psi::Array{T,2}
  sigma::Union(Array{T,2},Array{T,1}) # Innovation variance-covariance matrix
  grc::S                              # Number of eigenvalues greater than cutoff
  soln_type::String                   # "determinate", "indeterminate", or "explosive"
end
