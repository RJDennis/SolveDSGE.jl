function check_model_form(m)

  if typeof(m) <:Blanchard_Kahn_Form
    assessment = check_blanchard_kahn_form(m)
  elseif typeof(m) <:Klein_Form
    assessment = check_klein_form(m)
  elseif typeof(m) <:Structural_Form
    assessment = check_structural_form(m)
  elseif typeof(m) <:Sims_Form
    assessment = check_sims_form(m)
  elseif typeof(m) <:Gomme_Klein_Form
    assessment = check_gomme_klein_form(m)
  elseif typeof(m) <:Lombardo_Sutherland_Form
    assessment = check_lombardo_sutherland_form(m)
  end

  return assessment

end

function check_blanchard_kahn_form(m::Blanchard_Kahn_Form)

  assessment = "Pass"

  if m.nx < 1
    assessment = "Fail: Too few predetermined variables"
  elseif m.ny < 1
    assessment = "Fail: Too few non-predetermined variables"
  elseif m.nx+m.ny != size(m.a,1)
    assessment = "Fail: A is not conformable with nx and ny"
  elseif size(m.a,1) != size(m.a,2)
    assessment = "Fail: A is not square"
  elseif size(m.a,1) != size(m.c,1)
    assessment = "Fail: A and C must have the same number of rows"
  elseif size(m.c,2) > m.nx
    assessment = "Fail: There are more innovations than predetermined variables"
    elseif size(m.sigma,1) != size(m.sigma,2)
    assessment = "Fail: SIGMA must be square"
  elseif size(m.c,2) != size(m.sigma,1)
    assessment = "Fail: C and SIGMA are not conformable"
  end

  return assessment

end

function check_klein_form(m::Klein_Form)

  assessment = "Pass"

  if m.nx < 1
    assessment = "Fail: Too few predetermined variables"
  elseif m.ny < 1
    assessment = "Fail: Too few non-predetermined variables"
  elseif m.nx+m.ny != size(m.a,1)
    assessment = "Fail: A is not conformable with nx and ny"
  elseif size(m.a,1) != size(m.a,2)
    assessment = "Fail: A is not square"
  elseif size(m.b,1) != size(m.b,2)
    assessment = "Fail: B is not square"
  elseif size(m.a) != size(m.b)
    assessment = "Fail: A and B are not conformable"
  elseif size(m.a,1) != size(m.c,1)
    assessment = "Fail: A and C must have the same number of rows"
  elseif size(m.c,2) > m.nx
    assessment = "Fail: There are more innovations than predetermined variables"
    elseif size(m.sigma,1) != size(m.sigma,2)
    assessment = "Fail: SIGMA must be square"
  elseif size(m.c,2) != size(m.sigma,1)
    assessment = "Fail: C and SIGMA are not conformable"
  end

  return assessment

end

function check_structural_form(m::Structural_Form)

  assessment = "Pass"

  if size(m.a,1) != size(m.a,2)
    assessment = "Fail: A is not square"
  elseif size(m.a) != size(m.a1)
    assessment = "Fail: A and A1 are not conformable"
  elseif size(m.a) != size(m.b)
    assessment = "Fail: A and B are not conformable"
  elseif size(m.a,1) != size(m.c,1)
    assessment = "Fail: A and C must have the same number of rows"
  elseif size(m.sigma,1) != size(m.sigma,2)
    assessment = "Fail: SIGMA must be square"
  elseif size(m.c,2) != size(m.sigma,1)
    assessment = "Fail: C and SIGMA are not conformable"
  end

  return assessment

end

function check_sims_form(m::Sims_Form)

  assessment = "Pass"

  if size(m.gamma0,1) != size(m.gamma0,2)
    assessment = "Fail: GAMMA0 is not square"
  elseif size(m.gamma1,1) != size(m.gamma1,2)
    assessment = "Fail: GAMMA1 is not square"
  elseif size(m.gamma0) != size(m.gamma1)
    assessment = "Fail: GAMMA0 and GAMMA1 are not conformable"
  elseif size(m.gamma0,1) != size(m.c,1)
    assessment = "Fail: GAMMA0 and C must have the same number of rows"
  elseif size(m.c,2) != 1
    assessment = "Fail: C must be a vector"
  elseif size(m.gamma0,1) != size(m.psi,1)
    assessment = "Fail: GAMMA0 and PSI must have the same number of rows"
  elseif size(m.gamma0,1) > size(m.pi,1)
    assessment = "Fail: GAMMA0 and PI must have the same number of rows"
  elseif size(m.sigma,1) != size(m.sigma,2)
    assessment = "Fail: SIGMA must be square"
  elseif size(m.psi,2) != size(m.sigma,1)
    assessment = "Fail: PSI and SIGMA are not conformable"
  end

  return assessment

end

function check_gomme_klein_form(m::Gomme_Klein_Form)

  assessment = "Pass"

  if m.nx < 1
    assessment = "Fail: Too few predetermined variables"
  elseif m.ny < 1
    assessment = "Fail: Too few non-predetermined variables"
  elseif m.nx+m.ny != size(m.derivs,1)
    assessment = "Fail: Need first derivative for all model equations"
  elseif 2*size(m.derivs,1) != size(m.derivs,2)
    assessment = "Fail: The dimensions of the first-derivative matrix are incorrect"
  elseif size(m.derivs,2) != size(m.hessians,2)
    assessment = "Fail: The derivative matrix and the hessian matrix are not conformable"
  elseif size(m.hessians,1) != (m.nx+m.ny)*size(m.hessians,2)
    assessment = "Fail: The dimensions of the hessian matrix are incorrect"
  elseif size(m.eta,1) != m.nx
    assessment = "Fail: Eta has too few rows"
  elseif size(m.eta,2) != size(m.sigma,1)
    assessment = "Fail: Eta has too few columns"
  elseif size(m.sigma,1) != size(m.sigma,2)
    assessment = "Fail: SIGMA must be square"
  end

  return assessment

end

function check_lombardo_sutherland_form(m::Lombardo_Sutherland_Form)

  assessment = "Pass"

  if m.nx < 1
    assessment = "Fail: Too few predetermined variables"
  elseif m.ny < 1
    assessment = "Fail: Too few non-predetermined variables"
  elseif m.nx+m.ny != size(m.derivs,1)
    assessment = "Fail: Need first derivative for all model equations"
  elseif 2*size(m.derivs,1) != size(m.derivs,2)
    assessment = "Fail: The dimensions of the first-derivative matrix are incorrect"
  elseif size(m.derivs,2) != size(m.hessians,2)
    assessment = "Fail: The derivative matrix and the hessian matrix are not conformable"
  elseif size(m.hessians,1) != (m.nx+m.ny)*size(m.hessians,2)
    assessment = "Fail: The dimensions of the hessian matrix are incorrect"
  elseif size(m.eta,1) != m.nx
    assessment = "Fail: Eta has too few rows"
  elseif size(m.eta,2) != size(m.sigma,1)
    assessment = "Fail: Eta has too few columns"
  elseif size(m.sigma,1) != size(m.sigma,2)
    assessment = "Fail: SIGMA must be square"
  end

  return assessment

end
