function y = BEM_matvec (c, x, Gidx)
% Gidx = 1 (Gs), 2 (Gn)
  switch (c.cbem.use)
    case 2
      % hmmvp.
      y = hmmvp('mvp', c.cbem.id(Gidx), x);
    case 1
      % Old H-matrix.
      y = hsvd_mex(c.cbem.h(c.cbem.strt(Gidx) + c.cbem.thidx).id, x);
    case 0
      if (Gidx == 1) y = c.Gs*x;
      else           y = c.Gn*x; end
    otherwise
      error('c.cbem.use should be one of 0, 1, 2.');
  end
end
