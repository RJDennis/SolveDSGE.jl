function reorder_complex_schur{T<:AbstractFloat}(s::Array{T,2},q::Array{T,2})

  n = size(s,1)

  for ilast = n:-1:2

    mag = norm(s[ilast,ilast])
    icurrent = ilast

    for j = (ilast-1):-1:1

      if norm(s[j,j]) > mag
        mag = norm(s[j,j])
        icurrent = j
      end

    end

    for i = icurrent:(ilast-1)

      u = s[i:i+1,i+1]
      u[2] = u[2]-s[i,i]
      g = givens(u'[:])
      s[i:i+1,:] = g's[i:i+1,:]
      s[:,i:i+1] = s[:,i:i+1]*g
      q[i:i+1,:] = g'q[i:i+1,:]

    end

  end

  return s,q

end

function reorder_complex_schur{T<:AbstractFloat}(s::Array{Complex{T},2},q::Array{Complex{T},2})

  n = size(s,1)

  for ilast = n:-1:2

    mag = norm(s[ilast,ilast])
    icurrent = ilast

    for j = (ilast-1):-1:1

      if norm(s[j,j]) > mag
        mag = norm(s[j,j])
        icurrent = j
      end

    end

    for i = icurrent:(ilast-1)

      u = s[i:i+1,i+1]
      u[2] = u[2]-s[i,i]
      g = givens(u'[:])
      s[i:i+1,:] = g's[i:i+1,:]
      s[:,i:i+1] = s[:,i:i+1]*g
      q[i:i+1,:] = g'q[i:i+1,:]

    end

  end

  return s,q

end

function reorder_generalized_complex_schur{T<:AbstractFloat}(s::Array{T,2},t::Array{T,2},q::Array{T,2},z::Array{T,2})

  # This draws on draws on Kressner, D., (2006), "BLOCK ALGORITHMS FOR REORDERING STANDARD AND GENERALIZED SCHUR FORMS", LAPACK WORKING NOTE 171

  n = size(s,1)

  for ilast = n:-1:2

    if norm(t[ilast,ilast]) < 2*eps(T)*norm(t)
      mag = Inf
    else
      mag = norm(s[ilast,ilast]/t[ilast,ilast])
    end
    icurrent = ilast

    if mag != Inf
      for j = (ilast-1):-1:1

        if norm(t[j,j]) < 2*eps(T)*norm(t)
          mag = Inf
          icurrent = j
        elseif norm(s[j,j]/t[j,j]) > mag
          mag = norm(s[j,j]/t[j,j])
          icurrent = j
        end

      end
    end

    for i = icurrent:(ilast-1)

      s1 = s[i,i]
      s2 = s[i,i+1]
      s4 = s[i+1,i+1]
      t1 = t[i,i]
      t2 = t[i,i+1]
      t4 = t[i+1,i+1]

      if norm(t1) < 2*eps(T)*norm(t) && norm(s1) < 2*eps(T)*norm(s)
        error("Bz-A is not a regular matrix pencil: your model is not well posed")
      elseif norm(t4) < 2*eps(T)*norm(t) && norm(s4) < 2*eps(T)*norm(s)
        error("Bz-A is not a regular matrix pencil: your model is not well posed")
      elseif norm(t1) < 2*eps(T)*norm(t) && norm(t4) < 2*eps(T)*norm(s)
        # Both eigenvalues are infinite and do not need to be switched
        gq = eye(2)
        gz = eye(2)
      elseif norm(t1*s4-t4*s1) < 2*eps(T)*maximum([norm(s),norm(t)])
        # Both eigenvalues have equal modulus and do not need to be switched
        gq = eye(2)
        gz = eye(2)
      else
        h = [s1 -s4; t1 -t4]
        f = [s2; t2]
        p = inv(h)*f
        if norm(p[1]) < 2*eps(T)*norm(h) && norm(p[2]) < 2*eps(T)*norm(h)
          gq = [0.0 1.0; 1.0 0.0]
          gz = [0.0 1.0; 1.0 0.0]
        else
          u = [-p[2], complex(1.0)]
          (gq,r) = householder(u)
          u = [complex(1.0) p[1]]'
          (gz,r) = householder(u)
          gz = [gz[:,2] gz[:,1]]
        end
      end

      s[i:i+1,:] = gq's[i:i+1,:]
      t[i:i+1,:] = gq't[i:i+1,:]
      s[:,i:i+1] = s[:,i:i+1]*gz
      t[:,i:i+1] = t[:,i:i+1]*gz
      q[i:i+1,:] = gq'q[i:i+1,:]
      z[:,i:i+1] = z[:,i:i+1]*gz

    end

  end

  return s,t,q,z

end

function reorder_generalized_complex_schur{T<:AbstractFloat}(s::Array{Complex{T},2},t::Array{Complex{T},2},q::Array{Complex{T},2},z::Array{Complex{T},2})

  # This method draws on Kressner, D., (2006), "BLOCK ALGORITHMS FOR REORDERING STANDARD AND GENERALIZED SCHUR FORMS", LAPACK WORKING NOTE 171

  n = size(s,1)

  for ilast = n:-1:2

    if norm(t[ilast,ilast]) < 2*eps(T)*norm(t)
      mag = Inf
    else
      mag = norm(s[ilast,ilast]/t[ilast,ilast])
    end
    icurrent = ilast

    if mag != Inf
      for j = (ilast-1):-1:1

        if norm(t[j,j]) < 2*eps(T)*norm(t)
          mag = Inf
          icurrent = j
        elseif norm(s[j,j]/t[j,j]) > mag
          mag = norm(s[j,j]/t[j,j])
          icurrent = j
        end

      end
    end

    for i = icurrent:(ilast-1)

      s1 = s[i,i]
      s2 = s[i,i+1]
      s4 = s[i+1,i+1]
      t1 = t[i,i]
      t2 = t[i,i+1]
      t4 = t[i+1,i+1]

      if norm(t1) < 2*eps(T)*norm(t) && norm(s1) < 2*eps(T)*norm(s)
        error("Bz-A is not a regular matrix pencil: your model is not well posed")
      elseif norm(t4) < 2*eps(T)*norm(t) && norm(s4) < 2*eps(T)*norm(s)
        error("Bz-A is not a regular matrix pencil: your model is not well posed")
      elseif norm(t1) < 2*eps(T)*norm(t) && norm(t4) < 2*eps(T)*norm(s)
        # Both eigenvalues are infinite and do not need to be switched
        gq = eye(2)
        gz = eye(2)
      elseif norm(t1*s4-t4*s1) < 2*eps(T)*maximum([norm(s),norm(t)])
        # Both eigenvalues have equal modulus and do not need to be switched
        gq = eye(2)
        gz = eye(2)
      else
        h = [s1 -s4; t1 -t4]
        f = [s2; t2]
        p = inv(h)*f
        if norm(p[1]) < 2*eps(T)*norm(h) && norm(p[2]) < 2*eps(T)*norm(h)
          gq = [0.0 1.0; 1.0 0.0]
          gz = [0.0 1.0; 1.0 0.0]
        else
          u = [-p[2], complex(1.0)]
          (gq,r) = householder(u)
          u = [complex(1.0) p[1]]'
          (gz,r) = householder(u)
          gz = [gz[:,2] gz[:,1]]
        end
      end

      s[i:i+1,:] = gq's[i:i+1,:]
      t[i:i+1,:] = gq't[i:i+1,:]
      s[:,i:i+1] = s[:,i:i+1]*gz
      t[:,i:i+1] = t[:,i:i+1]*gz
      q[i:i+1,:] = gq'q[i:i+1,:]
      z[:,i:i+1] = z[:,i:i+1]*gz

    end

  end

  return s,t,q,z

end

#=

function reorder_complex_schur{T<:AbstractFloat}(s::Array{T,2},q::Array{T,2})

  # Alternative reordering using Householder transform.  This method draws on Kressner, D., (2006), "BLOCK ALGORITHMS FOR REORDERING
  # STANDARD AND GENERALIZED SCHUR FORMS", LAPACK WORKING NOTE 171

  n = size(s,1)

  for ilast = n:-1:2

    mag = norm(s[ilast,ilast])
    icurrent = ilast

    for j = (ilast-1):-1:1

      if norm(s[j,j]) > mag
        mag = norm(s[j,j])
        icurrent = j
      end

    end

    for i = icurrent:(ilast-1)

      s1 = s[i,i]
      s2 = s[i,i+1]
      s4 = s[i+1,i+1]

      if norm(s1-s4) < 2*eps(T)*norm(s)
        g = eye(2)
      else
        (g,r) = householder([-(s1-s4)^(-1)*s2,1.0])
      end

      s[i:i+1,:] = g's[i:i+1,:]
      s[:,i:i+1] = s[:,i:i+1]*g
      q[i:i+1,:] = g'q[i:i+1,:]

    end

  end

  return s,q

end

function reorder_complex_schur{T<:AbstractFloat}(s::Array{Complex{T},2},q::Array{Complex{T},2})

  # Alternative reordering using Householder transform.  This method draws on Kressner, D., (2006), "BLOCK ALGORITHMS FOR REORDERING
  # STANDARD AND GENERALIZED SCHUR FORMS", LAPACK WORKING NOTE 171

  n = size(s,1)

  for ilast = n:-1:2

    mag = norm(s[ilast,ilast])
    icurrent = ilast

    for j = (ilast-1):-1:1

      if norm(s[j,j]) > mag
        mag = norm(s[j,j])
        icurrent = j
      end

    end

    for i = icurrent:(ilast-1)

      s1 = s[i,i]
      s2 = s[i,i+1]
      s4 = s[i+1,i+1]

      if norm(s1-s4) < 2*eps(T)*norm(s)
        g = eye(2)
      else
        p = [-(s1-s4)^(-1)*s2,complex(1.0)]
        (g,r) = householder(p)
      end

      s[i:i+1,:] = g's[i:i+1,:]
      s[:,i:i+1] = s[:,i:i+1]*g
      q[i:i+1,:] = g'q[i:i+1,:]

    end

  end

  return s,q

end

=#
