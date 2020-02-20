function doubling(a::Array{T, 2}, b::Array{T, 2}, g::Array{T, 2}, tol::T, maxiters::S) where {T <: AbstractFloat, S <: Int}

  retcode = 0
  iters = 1
  len = Inf

  gn = similar(g)

  while len > tol && iters <= maxiters
    gn = g+a*g*b
    len = maximum(abs,gn-g)
    a = a*a
    b = b*b
    g = copy(gn)
    iters += 1
  end

  if iters > maxiters
    retcode = 1
  end

  return gn,retcode

end
