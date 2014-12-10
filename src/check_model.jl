function check_model_form(model)

  if typeof(model) <:Blanchard_Kahn_Form
    assessment = check_blanchard_kahn_form(model)
  elseif typeof(model) <:Klein_Form
    assessment = check_klein_form(model)
  elseif typeof(model) <:Binder_Pesaran_Form
    assessment = check_binder_pesaran_form(model)
  elseif typeof(model) <:Sims_Form
    assessment = check_sims_form(model)
  elseif typeof(model) <:Gomme_Klein_Form
    assessment = check_gomme_klein_form(model)
  elseif typeof(model) <:Lombardo_Sutherland_Form
    assessment = check_lombardo_sutherland_form(model)
  end

  return assessment

end

function check_blanchard_kahn_form(model::Blanchard_Kahn_Form)

  assessment = "Pass"

  if model.nx < 1
    assessment = "Fail: Too few predetermined variables"
  elseif model.ny < 1
    assessment = "Fail: Too few non-predetermined variables"
  elseif model.nx+model.ny != size(model.a,1)
    assessment = "Fail: A is not conformable with nx and ny"
  elseif size(model.a,1) != size(model.a,2)
    assessment = "Fail: A is not square"
  elseif size(model.a,1) != size(model.c,1)
    assessment = "Fail: A and C must have the same number of rows"
  elseif size(model.c,2) > model.nx
    assessment = "Fail: There are more innovations than predetermined variables"
  elseif size(model.sigma,1) != size(model.sigma,2)
    assessment = "Fail: SIGMA must be square"
  elseif size(model.c,2) != size(model.sigma,1)
    assessment = "Fail: C and SIGMA are not conformable"
  end

  return assessment

end

function check_klein_form(model::Klein_Form)

  assessment = "Pass"

  if model.nx < 1
    assessment = "Fail: Too few predetermined variables"
  elseif model.ny < 1
    assessment = "Fail: Too few non-predetermined variables"
  elseif model.nx+model.ny != size(model.a,1)
    assessment = "Fail: A is not conformable with nx and ny"
  elseif size(model.a,1) != size(model.a,2)
    assessment = "Fail: A is not square"
  elseif size(model.b,1) != size(model.b,2)
    assessment = "Fail: B is not square"
  elseif size(model.a) != size(model.b)
    assessment = "Fail: A and B are not conformable"
  elseif size(model.a,1) != size(model.c,1)
    assessment = "Fail: A and C must have the same number of rows"
  elseif size(model.c,2) > model.nx
    assessment = "Fail: There are more innovations than predetermined variables"
    elseif size(model.sigma,1) != size(model.sigma,2)
    assessment = "Fail: SIGMA must be square"
  elseif size(model.c,2) != size(model.sigma,1)
    assessment = "Fail: C and SIGMA are not conformable"
  end

  return assessment

end

function check_binder_pesaran_form(model::Binder_Pesaran_Form)

  assessment = "Pass"

  if size(model.a,1) != size(model.a,2)
    assessment = "Fail: A is not square"
  elseif size(model.a) != size(model.a1)
    assessment = "Fail: A and A1 are not conformable"
  elseif size(model.a) != size(model.b)
    assessment = "Fail: A and B are not conformable"
  elseif size(model.a,1) != size(model.c,1)
    assessment = "Fail: A and C must have the same number of rows"
  elseif size(model.sigma,1) != size(model.sigma,2)
    assessment = "Fail: SIGMA must be square"
  elseif size(model.c,2) != size(model.sigma,1)
    assessment = "Fail: C and SIGMA are not conformable"
  end

  return assessment

end

function check_sims_form(model::Sims_Form)

  assessment = "Pass"

  if size(model.gamma0,1) != size(model.gamma0,2)
    assessment = "Fail: GAMMA0 is not square"
  elseif size(model.gamma1,1) != size(model.gamma1,2)
    assessment = "Fail: GAMMA1 is not square"
  elseif size(model.gamma0) != size(model.gamma1)
    assessment = "Fail: GAMMA0 and GAMMA1 are not conformable"
  elseif size(model.gamma0,1) != size(model.c,1)
    assessment = "Fail: GAMMA0 and C must have the same number of rows"
  elseif size(model.c,2) != 1
    assessment = "Fail: C must be a vector"
  elseif size(model.gamma0,1) != size(model.psi,1)
    assessment = "Fail: GAMMA0 and PSI must have the same number of rows"
  elseif size(model.gamma0,1) > size(model.pi,1)
    assessment = "Fail: GAMMA0 and PI must have the same number of rows"
  elseif size(model.sigma,1) != size(model.sigma,2)
    assessment = "Fail: SIGMA must be square"
  elseif size(model.psi,2) != size(model.sigma,1)
    assessment = "Fail: PSI and SIGMA are not conformable"
  end

  return assessment

end

function check_gomme_klein_form(model::Gomme_Klein_Form)

  assessment = "Pass"

  if model.nx < 1
    assessment = "Fail: Too few predetermined variables"
  elseif model.ny < 1
    assessment = "Fail: Too few non-predetermined variables"
  elseif model.nx+model.ny != size(model.derivs,1)
    assessment = "Fail: Need first derivative for all model equations"
  elseif 2*size(model.derivs,1) != size(model.derivs,2)
    assessment = "Fail: The dimensions of the first-derivative matrix are incorrect"
  elseif size(model.derivs,2) != size(model.hessians,2)
    assessment = "Fail: The derivative matrix and the hessian matrix are not conformable"
  elseif size(model.hessians,1) != (model.nx+model.ny)*size(model.hessians,2)
    assessment = "Fail: The dimensions of the hessian matrix are incorrect"
  elseif size(model.eta,1) != model.nx
    assessment = "Fail: Eta has too few rows"
  elseif size(model.eta,2) != size(model.sigma,1)
    assessment = "Fail: Eta has too few columns"
  elseif size(model.sigma,1) != size(model.sigma,2)
    assessment = "Fail: SIGMA must be square"
  end

  return assessment

end

function check_lombardo_sutherland_form(model::Lombardo_Sutherland_Form)

  assessment = "Pass"

  if model.nx < 1
    assessment = "Fail: Too few predetermined variables"
  elseif model.ny < 1
    assessment = "Fail: Too few non-predetermined variables"
  elseif model.nx+model.ny != size(model.derivs,1)
    assessment = "Fail: Need first derivative for all model equations"
  elseif 2*size(model.derivs,1) != size(model.derivs,2)
    assessment = "Fail: The dimensions of the first-derivative matrix are incorrect"
  elseif size(model.derivs,2) != size(model.hessians,2)
    assessment = "Fail: The derivative matrix and the hessian matrix are not conformable"
  elseif size(model.hessians,1) != (model.nx+model.ny)*size(model.hessians,2)
    assessment = "Fail: The dimensions of the hessian matrix are incorrect"
  elseif size(model.eta,1) != model.nx
    assessment = "Fail: Eta has too few rows"
  elseif size(model.eta,2) != size(model.sigma,1)
    assessment = "Fail: Eta has too few columns"
  elseif size(model.sigma,1) != size(model.sigma,2)
    assessment = "Fail: SIGMA must be square"
  end

  return assessment

end
