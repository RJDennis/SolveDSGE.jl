function derivative2{T<:AbstractFloat}(f::Function,x::Array{T,1})

  n = length(f(x))
  m = length(x)

	e  = eps(T)^(1/3)*maxabs([x;one(T)])
  dh = eye(m)*e
	deriv = Array(T,n,m)

    for i = 1:m
      deriv[:,i] = (-f(x+2*dh[:,i])+8*f(x+dh[:,i])-8*f(x-dh[:,i])+f(x-2*dh[:,i]))/(12*e)
    end

    return deriv

end
