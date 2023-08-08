function ok = makemex(p)
  if(nargin < 1) p = '.'; end
  flags = '';
  flags = ['CFLAGS="\$CFLAGS -fexceptions  -std=c99 -fopenmp -Wall" ',...
          'LDFLAGS="\$LDFLAGS -fopenmp"'];
  flags2 = ['-lmwblas -lmwlapack'];
  mc = sprintf('mex -R2017b -O -outdir %s %s',p,flags);
  fns = {'hsvd_mex.c'};
  try
    for(i = 1:length(fns))
      eval(sprintf('%s %s/%s %s/utils.c %s',mc,p,fns{i},p,flags2));
    end
    ok=1;
  catch
    ok=0;
  end  
end
