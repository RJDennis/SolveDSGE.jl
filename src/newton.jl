function newton(f::Function,x::Array{Float64,1},tol::Float64,maxiters::Int64)

  n = length(x)
  x_new = similar(x)

  iter = 0
  len = Inf
  while len > tol

    g = derivative(f,x)
  	update = g\f(x)
  	for i = 1:n
  	  x_new[i] = x[i] - update[i]
  	end

    len = maximum(abs,x_new-x)
    x = copy(x_new)

    iter += 1
    if iter >= maxiters
      break
    end

  end

  return x_new, f(x_new), iter

end
