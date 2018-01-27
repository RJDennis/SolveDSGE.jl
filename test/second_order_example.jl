# This code demonstrates the use of SolveDSGE for solving rational expectations models
# to second-order accuracy.  The demonstration models are the first and second examples
# from Schmitt-Grohe and Uribe (2004)

#-------------------------------------------------------------------------------------------#

# Bring the SolveDSGE module into global scope

using SolveDSGE

# Bring the NLsolve module into global scope

# Example one from Schmitt-Grohe and Uribe (2004)

# Specify the model in terms of a function whose fix-point in the model's nonstochastic steady state

function model_eqns(x)

  f = zeros(length(x))

  f[1] = x[4] - rho*x[1]
  f[2] = (1.0-delta)*exp(x[2]) + exp(x[1])*exp(x[2])^alpha - exp(x[3]) - exp(x[5])
  f[3] = exp(x[3])^(-gama) - betta*exp(x[6])^(-gama)*(1.0-delta+alpha*exp(x[4])*exp(x[5])^(alpha-1.0))
  f[4] = x[4] - x[1]
  f[5] = x[5] - x[2]
  f[6] = x[6] - x[3]

#  println(x)

  return f

end

# Specify the model in terms of a function to compute its first derivatives
# The function above (model_eqns) can also be used for this purpose

function model_eqns_2(x)

  f = zeros(3)

  f[1] = x[4] - rho*x[1]
  f[2] = (1.0-delta)*exp(x[2]) + exp(x[1])*exp(x[2])^alpha - exp(x[3]) - exp(x[5])
  f[3] = exp(x[3])^(-gama) - betta*exp(x[6])^(-gama)*(1.0-delta+alpha*exp(x[4])*exp(x[5])^(alpha-1.0))

  return f

end

# Specify the three model equations separately in order to compute the second derivatives

# Technology shock equation

function model_eqn_1(x)

  f = [0.0]

  f[1] = x[4] - rho*x[1]

  return f

end

# Capital accumulation equation

function model_eqn_2(x)

  f = [0.0]

  f[1] = (1.0-delta)*exp(x[2]) + exp(x[1])*exp(x[2])^alpha - exp(x[3]) - exp(x[5])

  return f

end

# Consumption Euler equation

function model_eqn_3(x)

  f = [0.0]

  f[1] = exp(x[3])^(-gama) - betta*exp(x[6])^(-gama)*(1.0-delta+alpha*exp(x[4])*exp(x[5])^(alpha-1.0))

  return f

end

# Assign values to model parameters

betta  = 0.95
alpha = 0.3
gama  = 2.0
delta = 1.0
rho   = 0.0
etta   = [1.0, 0.0]
sigma = [1.0]

# Initial guess for computing steady state

x = [0.01, -1.8, -0.9]
x = [x; x]

# Solve for the nonstochastic steady state

(z, f_zero, iters) = newton(model_eqns, x, 1e-14, 1000)

ss = copy(z)

# Compute the model's first derivatives at the nonstochastic steady state

deriv1 = derivative(model_eqns_2,ss)

# Compute the model's second derivatives at the nonstochastic steady state

deriv21 = hessian(model_eqn_1,ss)
deriv22 = hessian(model_eqn_2,ss)
deriv23 = hessian(model_eqn_3,ss)

deriv2 = [deriv21; deriv22; deriv23]

nx = 2
ny = 1
cutoff = 1.0

# Put the model into Gomme_Klein_Form type

m_gk = Gomme_Klein_Form(nx,ny,deriv1,deriv2,etta,sigma)

# Do some basic checking to determine whether the model matrices are conformable

assessment = check_model_form(m_gk)
println(assessment)

# Solve the model using Gomme and Klein (2011)

soln_gk = solve_re(m_gk,cutoff)
responses_gk = impulses(soln_gk,5,1)

# Put the model into Lombardo_Sutherland_Form type

m_ls = Lombardo_Sutherland_Form(nx,ny,deriv1,deriv2,etta,sigma)

# Solve the model using Lombardo and Sutherland (2007)

soln_ls = solve_re(m_ls,cutoff)
responses_ls_1 = impulses(soln_ls,5,1)

# Convert Lombardo_Sutherland_Soln to Gomme_Klein_Soln

new_soln_ls = convert_second_order(soln_ls)
responses_ls_2 = impulses(new_soln_ls,5,1)

# Check to see that the two solutions are the same

println(maximum(abs,new_soln_ls.hxx-soln_gk.hxx))
println(maximum(abs,new_soln_ls.gxx-soln_gk.gxx))

println(maximum(abs,responses_gk-responses_ls_1))

# Example two from Schmitt-Grohe and Uribe (2004)

# Specify the model in terms of a function whose fix-point in the model's nonstochastic steady state

function model2_eqns(x)

  f = zeros(length(x))

  f[1]  = x[7] - rho*x[1]
  f[2]  = x[8] - rho*x[2]
  f[3]  = exp(x[5]) + exp(x[6]) + exp(x[9]) - (1.0-delta)*exp(x[3]) + exp(x[10]) - (1.0-delta)*exp(x[4]) - exp(x[1])*exp(x[3])^alpha - exp(x[2])*exp(x[4])^alpha
  f[4]  = exp(x[5]) - exp(x[6])
  f[5]  = exp(x[5])^(-gama) - betta*exp(x[11])^(-gama)*(1.0-delta+alpha*exp(x[7])*exp(x[9])^(alpha-1.0))
  f[6]  = exp(x[6])^(-gama) - betta*exp(x[12])^(-gama)*(1.0-delta+alpha*exp(x[8])*exp(x[10])^(alpha-1.0))
  f[7]  = x[7] - x[1]
  f[8]  = x[8] - x[2]
  f[9]  = x[9] - x[3]
  f[10] = x[10] - x[4]
  f[11] = x[11] - x[5]
  f[12] = x[12] - x[6]

  return f

end

# Specify the model in terms of a function to compute its first derivatives
# The function above (model2_eqns) can also be used for this purpose

function model2_eqns_2(x)

  f = zeros(length(x))

  f[1]  = x[7] - rho*x[1]
  f[2]  = x[8] - rho*x[2]
  f[3]  = exp(x[5]) + exp(x[6]) + exp(x[9]) - (1.0-delta)*exp(x[3]) + exp(x[10]) - (1.0-delta)*exp(x[4]) - exp(x[1])*exp(x[3])^alpha - exp(x[2])*exp(x[4])^alpha
  f[4]  = exp(x[5]) - exp(x[6])
  f[5]  = exp(x[5])^(-gama) - betta*exp(x[11])^(-gama)*(1.0-delta+alpha*exp(x[7])*exp(x[9])^(alpha-1.0))
  f[6]  = exp(x[6])^(-gama) - betta*exp(x[12])^(-gama)*(1.0-delta+alpha*exp(x[8])*exp(x[10])^(alpha-1.0))

  return f

end

# Specify the three model equations separately in order to compute the second derivatives

# First technology shock equation

function model2_eqn_1(x)

  f = [0.0]

  f[1]  = x[7] - rho*x[1]

  return f

end

# Second technology shock equation

function model2_eqn_2(x)

  f = [0.0]

  f[1]  = x[8] - rho*x[2]

  return f

end

# Capital accumulation equation

function model2_eqn_3(x)

  f = [0.0]

  f[1]  = exp(x[5]) + exp(x[6]) + exp(x[9]) - (1.0-delta)*exp(x[3]) + exp(x[10]) - (1.0-delta)*exp(x[4]) - exp(x[1])*exp(x[3])^alpha - exp(x[2])*exp(x[4])^alpha

  return f

end

# Simple way of specifying second capital accumulation equation

function model2_eqn_4(x)

  f = [0.0]

  f[1]  = exp(x[5]) - exp(x[6])

  return f

end

# Euler equation for first consumption good

function model2_eqn_5(x)

  f = [0.0]

  f[1]  = exp(x[5])^(-gama) - betta*exp(x[11])^(-gama)*(1.0-delta+alpha*exp(x[7])*exp(x[9])^(alpha-1.0))

  return f

end

# Euler equation for second consumption good

function model2_eqn_6(x)

  f = [0.0]

  f[1]  = exp(x[6])^(-gama) - betta*exp(x[12])^(-gama)*(1.0-delta+alpha*exp(x[8])*exp(x[10])^(alpha-1.0))

  return f

end

# Assign values to model parameters

betta = 0.95
alpha = 0.3
delta = 0.1
gama  = 2.0
rho   = 0.0
etta   = [1.0 0.0; 0.0 1.0; 0.0 0.0; 0.0 0.0]
sigma = [1.0 0.0; 0.0 1.0]

# Initial guess for computing steady state

x = [0.1, 0.1, 1.0, 1.0, 0.1, 0.1]
x = [x; x]

# Solve for the nonstochastic steady state
(z, f_zero, iters) = newton(model2_eqns, x, 1e-14, 1000)

ss = copy(z)

# Compute model's first derivatives at steady state

deriv1 = derivative(model2_eqns_2,ss)

# Compute model's second derivatives at steady state

deriv21 = hessian(model2_eqn_1,ss)
deriv22 = hessian(model2_eqn_2,ss)
deriv23 = hessian(model2_eqn_3,ss)
deriv24 = hessian(model2_eqn_4,ss)
deriv25 = hessian(model2_eqn_5,ss)
deriv26 = hessian(model2_eqn_6,ss)

deriv2 = [deriv21; deriv22; deriv23; deriv24; deriv25; deriv26]

nx = 4
ny = 2
cutoff = 1.0

# Put the model into Gomme_Klein_Form type

m_gk = Gomme_Klein_Form(nx,ny,deriv1,deriv2,etta,sigma)

# Do some basic checking to determine whether the model matrices are conformable

assessment = check_model_form(m_gk)
println(assessment)

# Solve the model using Gomme and Klein (2011)

soln_gk = solve_re(m_gk,cutoff)
responses_gk = impulses(soln_gk,5,1)

# Put the model into Lombard_Sutherland_Form type

m_ls = Lombardo_Sutherland_Form(nx,ny,deriv1,deriv2,etta,sigma)

# Solve the model using Lombardo and Sutherland (2007)

soln_ls = solve_re(m_ls,cutoff)
responses_ls_1 = impulses(soln_ls,5,1)

# Convert Lombardo_Sutherland_Soln to Gomme_Klein_Soln

new_soln_ls = convert_second_order(soln_ls)
responses_ls_2 = impulses(new_soln_ls,5,1)

# Check to see that the two solutions are the same

println(maximum(abs,new_soln_ls.hxx-soln_gk.hxx))
println(maximum(abs,new_soln_ls.gxx-soln_gk.gxx))

println(maximum(abs,responses_gk-responses_ls_1))
