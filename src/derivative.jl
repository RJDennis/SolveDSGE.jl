function derivative{T<:AbstractFloat}(f::Function,x::Array{T,1})

    n = length(f(x))
    m = length(x)

	e  = eps(T)^(1/3)*Base.maxabs([x,one(T)])
    dh = eye(m)*e
	deriv = Array(T,n,m)

    for i = 1:m
	    f1 = f(x+2*dh[:,i])
	    f2 = f(x+dh[:,i])
	    f3 = f(x-dh[:,i])
	    f4 = f(x-2*dh[:,i])
        deriv[:,i] = (-f1+8*f2-8*f3+f4)/(12*e)
    end

    return deriv

end
