function convert_second_order(model::M) where M <: Lombardo_Sutherland_Soln

  ssh       = copy(model.ssh)
  hx        = copy(model.hx)
  hxx       = copy(model.hxx)
  ssg       = copy(model.ssg)
  gx        = copy(model.gx)
  gxx       = copy(model.gxx)
  eta       = copy(model.eta)
  sigma     = copy(model.sigma)
  grc       = copy(model.grc)
  soln_type = model.soln_type

  (n1,n2) = size(hx)
  (n3,n4) = size(gx)

  new_hxx = reshape(vech_to_vec(n2)*hxx[1,:],n2,n2)
  for i = 2:n1
    new_hxx = [new_hxx; reshape(vech_to_vec(n2)*hxx[i,:],n2,n2)]
  end

  (n5,n6) = size(new_hxx)

  j = 0
  for i = 1:n5
    j += 1
    new_hxx[i,j] = 2*new_hxx[i,j]
    if j == n6
      j = 0
    end
  end

  new_hxx = new_hxx/2

  new_gxx = reshape(vech_to_vec(n4)*gxx[1,:],n4,n4)
  for i = 2:n3
    new_gxx = [new_gxx; reshape(vech_to_vec(n4)*gxx[i,:],n4,n4)]
  end

  (n7,n8) = size(new_gxx)

  j = 0
  for i = 1:n7
    j += 1
    new_gxx[i,j] = 2*new_gxx[i,j]
    if j == n8
      j = 0
    end
  end

  new_gxx = new_gxx/2

  soln = Gomme_Klein_Soln(ssh,hx,new_hxx,ssg,gx,new_gxx,eta,sigma,grc,soln_type)

  return soln

end
