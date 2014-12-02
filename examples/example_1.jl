# This code demonstrates the use of SolveDSGE for solving linear rational expectations models.
# The demonstration model is a simple forward-looking new Keynesian model with monetary policy
# conducted using a Taylor rule.

# The model is as follows:

# pi_{t} = beta*E_{t}pi_{t+1} + kappa*y_{t} + u_{t}
#  y_{t} = E_{t}y_{t+1} - mu*(R_{t}-E_{t}pi_{t+1} + v_{t}
#  R_{t} = phi_pi*E_{t}pi_{t+1} + phi_y*y_{t}

#-------------------------------------------------------------------------------------------#

# Bring the SolveDSGE module into global scope

using SolveDSGE

# Assign values to model parameters

const beta = 0.99
const mu = 0.5
const phi_pi = 1.5
const rho_u = 0.8
const rho_v = 0.8
const kappa = 0.1
const phi_y = 0.5
const sigma = [1.0 0.0; 0.0 1.0]

const cutoff = 1.0
const tol    = 1e-8

# Solve the model using Blanchard and Kahn (1980)

# Construct the model matrices

b = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 beta 0.0; 0.0 0.0 mu*(1.0-phi_pi) 1.0]
a = [rho_u 0.0 0.0 0.0; 0.0 rho_v 0.0 0.0; -1.0 0.0 1.0 -kappa; 0.0 -1.0 0.0 (1.0+mu*phi_y)]
c = [1.0 0.0; 0.0 1.0; 0.0 0.0; 0.0 0.0]

c = b\c
a = b\a

nx = 2
ny = 2

# Put the model into Blanchard_Kahn_Form type

m_bk = Blanchard_Kahn_Form(nx,ny,a,c,sigma)

# Do some basic checking to determine whether the model matrices are conformable

assessment = check_model_form(m_bk)
println(assessment)

# Solve the model by solving an eigenvalue problem

soln_bk_1 = solve_re(m_bk,cutoff)

# Solve the model by solving an algebraic Riccati equation

soln_bk_2 = solve_re(m_bk,cutoff,tol)

# Solve the model using Klein (2000)

# Construct the model matrices

b = [1.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 beta 0.0 0.0; 0.0 0.0 mu 1.0 0.0; 0.0 0.0 phi_pi 0.0 0.0]
a = [rho_u 0.0 0.0 0.0 0.0; 0.0 rho_v 0.0 0.0 0.0; -1.0 0.0 1.0 -kappa 0.0; 0.0 -1.0 0.0 1.0 mu; 0.0 0.0 0.0 -phi_y 1.0]
c = [1.0 0.0; 0.0 1.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]

nx = 2
ny = 3

# Put the model into Klein_Form type

m_k = Klein_Form(nx,ny,a,b,c,sigma)

# Do some basic checking to determine whether the model matrices are conformable

assessment = check_model_form(m_k)
println(assessment)

# Solve the model by solving a generalized eigenvalue problem

soln_k = solve_re(m_k,cutoff)

# Solve the model using Sims (2001)

# Construct the model matrices

gamma0 = [1.0 0.0 0.0 0.0 0.0 0.0 0.0;
          0.0 1.0 0.0 0.0 0.0 0.0 0.0;
          0.0 0.0 1.0 0.0 0.0 0.0 0.0;
          0.0 0.0 0.0 1.0 0.0 0.0 0.0;
          0.0 0.0 0.0 -phi_y 1.0 -phi_pi 0.0;
          -1.0 0.0 1.0 -kappa 0.0 -beta 0.0;
          0.0 -1.0 0.0 1.0 mu -mu -1.0]

gamma1 = [rho_v 0.0 0.0 0.0 0.0 0.0 0.0;
          0.0 rho_u 0.0 0.0 0.0 0.0 0.0;
          0.0 0.0 0.0 0.0 0.0 1.0 0.0;
          0.0 0.0 0.0 0.0 0.0 0.0 1.0;
          0.0 0.0 0.0 0.0 0.0 0.0 0.0;
          0.0 0.0 0.0 0.0 0.0 0.0 0.0;
          0.0 0.0 0.0 0.0 0.0 0.0 0.0]

c = zeros(7)

psi = [1.0 0.0;
       0.0 1.0;
       0.0 0.0;
       0.0 0.0;
       0.0 0.0;
       0.0 0.0;
       0.0 0.0]

pi = [0.0 0.0;
      0.0 0.0
      1.0 0.0
      0.0 1.0;
      0.0 0.0;
      0.0 0.0;
      0.0 0.0]

# Put the model into Sims_Form type

m_sims = Sims_Form(gamma0,gamma1,c,psi,pi,sigma)

# Do some basic checking to determine whether the model matrices are conformable

assessment = check_model_form(m_sims)
println(assessment)

# Solve the model by solving a generalized eigenvalue problem

soln_sims = solve_re(m_sims,cutoff)

# Solve the model in structural form

# Construct the model matrices

a0 = [1.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0; -1.0 0.0 1.0 -kappa 0.0; 0.0 -1.0 0.0 1.0 mu; 0.0 0.0 0.0 -phi_y 1.0]
a1 = [rho_u 0.0 0.0 0.0 0.0; 0.0 rho_v 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0]
a2 = [0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 beta 0.0 0.0; 0.0 0.0 mu 1.0 0.0; 0.0 0.0 phi_pi 0.0 0.0]
a4 = [1.0 0.0; 0.0 1.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]

# Put the model in Structural_Form type

m_s = Structural_Form(a0,a1,a2,a4,sigma)

# Do some basic checking to determine whether the model matrices are conformable

assessment = check_model_form(m_s)
println(assessment)

# Solve the by by solving a generalized eigenvalue problem

soln_s_1 = solve_re(m_s,cutoff)

# Solve the model using Binder and Pesaran's (1995) brute force method

soln_s_2 = solve_re(m_s,cutoff,tol)
