function hessian{T<:AbstractFloat}(f::Function,x::Array{T,1})

    m = length(x)

	e = eps(T)^(1/4)*maxabs([x;one(T)])
	dh = eye(m)*e
    hess = Array(T,m,m)

    for i = 1:m
	    f1 = f(x+2*dh[:,i])
		f2 = f(x+dh[:,i])
		f3 = f(x)
		f4 = f(x-dh[:,i])
		f5 = f(x-2*dh[:,i])
		hess[i:i,i:i] = (-f1+16*f2-30*f3+16*f4-f5)/(12*e^2)
    end

    for i = 1:m-1
        for j = i+1:m
		    f1 = f(x+dh[:,i]+dh[:,j])
			f2 = f(x-dh[:,i]-dh[:,j])
			f3 = f(x-dh[:,i]+dh[:,j])
			f4 = f(x+dh[:,i]-dh[:,j])
			hess[i:i,j:j] = (f1+f2-f3-f4)/(4*e^2)
            hess[j,i] = hess[i,j]
        end
    end

    return hess

end
