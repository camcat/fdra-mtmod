function [f f_psi f_gamma] = Friction(c,psi,gamma)
  deriv = nargout > 1;
  g = 0.5*exp(psi + (c.mu_0 + c.b.*gamma)./c.a);
  f = c.a.*asinh(g);
  if(deriv)
    gasinh_g = g./sqrt(g.^2 + 1);
    f_psi    = c.a.*gasinh_g;
    f_gamma  = gasinh_g.*c.b;
  end
