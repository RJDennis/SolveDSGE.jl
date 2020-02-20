using Pkg
Pkg.add("https://github.com/RJDennis/ChebyshevApprox.jl")
Pkg.add("https://github.com/RJDennis/SmolyakApprox.jl")
Pkg.add("https://github.com/RJDennis/PiecewiseLinearApprox.jl")

include("old/linear_example.jl")
include("old/second_order_example.jl")
