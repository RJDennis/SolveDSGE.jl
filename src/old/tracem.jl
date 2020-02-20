function tracem(x::Array{T, 2}) where T <: AbstractFloat

  # We require the number of rows to be greater than the number of columns, so that m is greater than one.

  trans = false

  if size(x,1) < size(x,2)
    x = x'
    trans = true
  end

  n = size(x,2)
  m = round(Int,size(x,1)/n)

  y = zeros(m,1)  # We want this to be a 2-d array for subsequent matrix multiplication

  for i = 1:m
    y[i,1] = tr(x[(n*(i-1)+1):i*n,1:n])
  end

  if trans == true
    y = y'
  end

  return y

end
