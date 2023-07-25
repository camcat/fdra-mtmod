function phi_dot = EvolveDlt (c, dlt, phi, psi, v, es)
%function phi_dot = EvolveDlt (c, dlt, phi, psi, v, es)
  % -theta_dot / theta
  chi = psi + dlt;
  switch (c.evolution)
   case 1 % aging law
     phi_dot = (1 - exp(-chi)).*v./c.d_c;
   case 2 % slip law
     phi_dot = chi.*v./c.d_c;
  end
end
