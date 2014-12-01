function vech_to_vec{S<:Int}(n::S)

  a = zeros(Int,n,n)
  tracker = 1

  for r = 1:n

    for c = 1:n

      if c >= r
        a[r,c] = tracker
        tracker += 1
      else
        a[r,c] = a[c,r]
      end

    end

  end

  gi = vec(a)

  x = zeros(Int,n*n,int(n*(n+1)/2))

  for r = 1:size(x,1)
    x[r,gi[r]] = 1
  end

  return x

end
