function vec_to_vech{S<:Int}(n::S)

  x = zeros(Int,round(Int,(n*(n+1)/2)),n*n)

  c = 1
  r = 1
  counter = 1

  for i = 1:n*n

    if r >= c
      x[counter,i] = 1
      counter += 1
    end

    r += 1
    if r > n
      c += 1
      r = 1
    end

  end

  return x

end
