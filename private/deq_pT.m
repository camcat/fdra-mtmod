function [dydt] = deq_pT(t,y,c)
% ODE function for use with ode23.

  global gP;
  dydt = zeros(size(y));
  
  slip    = y(1:c.Ncell);
  thetavd = y(c.Ncell+1:2*c.Ncell); % v_0 theta / d_c
  gamma   = log(thetavd);
  v       = c.v0*y(2*c.Ncell+1:3*c.Ncell);
  gP.psi  = log( y(2*c.Ncell+1:3*c.Ncell));
  
  % Stress
  tau = [];
  s = [slip(:); c.v_creep*(gP.t_g + t)];
  tau = c.tau0 + BEM_matvec(c,s,1); 
  gP.bes = c.bes; 
  
  phi_dot = EvolveDlt(c, gamma, [], gP.psi, v, gP.bes);

  s_eff = gP.bes;
 
  [mu mu_psi mu_gamma] = Friction(c,gP.psi,gamma);
  
  v1 = [v(:); c.v_creep];
  tau_dot = BEM_matvec(c,v1,1);
  c2 = 1./(mu_psi.*s_eff + c.eta.*v);
  q = tau_dot + s_eff.*mu_gamma.*phi_dot;
  gP.dvdt = v.*c2.*q;

  dydt = [v; -thetavd.*phi_dot; gP.dvdt/c.v0];
end

