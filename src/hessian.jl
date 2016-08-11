function hessian{T<:AbstractFloat}(f::Function,x::Array{T,1})

  m = length(x)

	e = eps(T)^(1/4)*maxabs([x;one(T)])
	dh = eye(m)*e
  hess = Array(T,m,m)

  for i = 1:m
		hess[i,i] = (-f(x+2*dh[:,i])+16*f(x+dh[:,i])-30*f(x)+16*f(x-dh[:,i])-f(x-2*dh[:,i]))/(12*e^2)
  end

    for i = 1:m-1
      for j = i+1:m
			  hess[i,j] = (f(x+dh[:,i]+dh[:,j])+f(x-dh[:,i]-dh[:,j])-f(x-dh[:,i]+dh[:,j])-f(x+dh[:,i]-dh[:,j]))/(4*e^2)
        hess[j,i] = hess[i,j]
      end
    end

    return hess

end
