import Base.givens

function givens(u::Array{T,1}) where T<:AbstractFloat

  if length(u) != 2
    error("givens: u should contain just two elements")
  end

  x = u[1]
  y = u[2]

  if y == 0.0
    c = 1.0
    s = 0.0
  elseif norm(y) >= norm(x)
    tau = -x/y
    s = 1.0/(sqrt(1.0+tau'tau))
    c = s*tau
  else
    tau = -y/x
    c = 1.0/(sqrt(1.0+tau'tau))
    s = c*tau
  end

  g = [c' s; -s' c]

  return g

end

function givens(u::Array{Complex{T},1}) where T<:AbstractFloat

  if length(u) != 2
    error("givens: u should contain just two elements")
  end

  x = u[1]
  y = u[2]

  if y == 0.0
    c = 1.0
    s = 0.0
  elseif norm(y) >= norm(x)
    tau = -x/y
    s = 1.0/(sqrt(1.0+tau'tau))
    c = s*tau
  else
    tau = -y/x
    c = 1.0/(sqrt(1.0+tau'tau))
    s = c*tau
  end

  g = [c' s; -s' c]

  return g

end
