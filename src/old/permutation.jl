function permutation(n1::S, n2::S) where S <: Int

  nn = n1*n2
  p = zeros(S,nn,nn)

  i = 1
  j = 1
  for r = 1:nn
    if j > n2
      j = 1
      i += 1
    end
    c = i+n1*(j-1)
    p[r,c] = 1

    j += 1
  end

  return p

end
