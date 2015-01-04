type State_Space_Objective{T<:FloatingPoint}
  beta::T
  q::Array{T,2}
  u::Array{T,2}
  r::Array{T,2}
end

type State_Space_Form{T<:FloatingPoint,S<:Integer}
  # E_tY[t+1] = A*Y[t] + B*U[t] + C*V[t+1]
  nx::S             # Number of predetermined variables
  ny::S             # Number of nonpredetermined variables
  a::Array{T,2}     # Companion matrix
  b::Array{T,2}     # Policy loadings
  c::Array{T,2}     # Innovation loading matrix
  sigma::Array{T,2} # Innovation variance-covariance matrix
end

type Generalized_State_Space_Form{T<:FloatingPoint,S<:Integer}
  # E_tY[t+1] = A*Y[t] + B*U[t] + C*V[t+1]
  nx::S             # Number of predetermined variables
  ny::S             # Number of nonpredetermined variables
  a0::Array{T,2}    # Companion matrix
  a::Array{T,2}     # Companion matrix
  b::Array{T,2}     # Policy loadings
  c::Array{T,2}     # Innovation loading matrix
  sigma::Array{T,2} # Innovation variance-covariance matrix
end

type State_Space_Soln{T<:FloatingPoint}
  p::Array{T,2}     # Transition matrix for predetermined variables
  k::Array{T,2}     # Innovation loading matrix
  h::Array{T,2}
  f::Array{T,2}     # Decision rule matrix linking nonpredetermined variables to predetermined variables
  v::Array{T,2}
  sigma::Array{T,2} # Innovation variance-covariance matrix
  converged::Bool
end

type Structural_Objective{T<:FloatingPoint}
  beta::T
  q::Array{T,2}
  r::Array{T,2}
end

type Structural_Form{T<:FloatingPoint}
  # E_tY[t+1] = A*Y[t] + B*U[t] + C*V[t+1]
  a0::Array{T,2}     # Companion matrix
  a1::Array{T,2}     # Policy loadings
  a2::Array{T,2}     # Innovation loading matrix
  a3::Array{T,2}     # Policy loadings
  a5::Array{T,2}     # Innovation loading matrix
  sigma::Array{T,2}  # Innovation variance-covariance matrix
end

type Generalized_Structural_Form{T<:FloatingPoint}
  # E_tY[t+1] = A*Y[t] + B*U[t] + C*V[t+1]
  a0::Array{T,2}     # Companion matrix
  a1::Array{T,2}     # Policy loadings
  a2::Array{T,2}     # Innovation loading matrix
  a3::Array{T,2}     # Policy loadings
  a4::Array{T,2}     # Policy loadings
  a5::Array{T,2}     # Innovation loading matrix
  sigma::Array{T,2}  # Innovation variance-covariance matrix
end

type Structural_Soln{T<:FloatingPoint}
  h::Array{T,2}     # Transition matrix for predetermined variables
  g::Array{T,2}     # Innovation loading matrix
  v::Array{T,2}
  sigma::Array{T,2} # Innovation variance-covariance matrix
  converged::Bool
end
