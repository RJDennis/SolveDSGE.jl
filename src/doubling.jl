function doubling{T<:AbstractFloat,S<:Int}(a::Array{T,2},b::Array{T,2},g::Array{T,2},tol::T,maxiters::S)

  retcode = 0
  j = 1
  len = Inf

  gn = copy(g)

  while len > tol && j <= maxiters
    an = a*a
    bn = b*b
    gn = g+a*g*b
    len = maxabs(gn-g)
    a = an
    b = bn
    g = copy(gn)
    j += 1
  end

  if j > maxiters
    retcode = 1
  end

  return gn,retcode

end
