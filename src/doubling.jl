function doubling{T<:AbstractFloat,S<:Int}(a::Array{T,2},b::Array{T,2},g::Array{T,2},tol::T,maxiters::S)

  retcode = 0
  j = 1
  len = Inf

  gn = zeros(g)

  while len > tol && j <= maxiters
    gn = g+a*g*b
    len = maxabs(gn-g)
    a = a*a
    b = b*b
    g = copy(gn)
    j += 1
  end

  if j > maxiters
    retcode = 1
  end

  return gn,retcode

end
