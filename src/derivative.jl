function derivative(f::Function, x::Array{T, 1}) where T <: AbstractFloat

  n = length(f(x))
  m = length(x)

	e  = eps(T)^(1/3)*maximum(abs,[x;one(T)])
  dh = Matrix{Float64}(I, m, m)*e
	deriv = Array{T}(undef,n,m)

  for i = 1:m
    deriv[:,i] = (-f(x+2*dh[:,i])+8*f(x+dh[:,i])-8*f(x-dh[:,i])+f(x-2*dh[:,i]))/(12*e)
  end

  return deriv

end
