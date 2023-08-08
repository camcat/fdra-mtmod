function varargout = Util(fn,varargin)
% In all these functions, you may omit '.dat' at the end of file names.
%
% s = qload(fn,skip,doStride)
%     Load a simulation using fewer keystrokes. The arguments are the same as
%   those for Read.
%     s.lo contains the load-mat struct. s.c contains FDRA-formatted physical
%   constants. s has several fields:
%     x: cell-centered along-fault coordinates
%     slip,v: slip and slip rate
%     psi: log(v/v0)
%     gam: log(v0 theta/d_c)
%     chi: psi + gam
%     (b)es: (background) effective stress.
%     mr: average slip speed over the fault
%     mu: friction coefficient
%     tau: transform of slip (and so includes radiation damping; compute
%       mu.*es to exclude damping)
%     p/Tprof: p and T profiles corresponding to s.lo.profs.
%
% [t psi chi dlt slip] =
%   Read(datFn,[stride],[doStride])
%   Read a .dat file. Optionally provide a stride number on the
%   columns. (Default is 1.) To read only, say, t and psi, omit the remaining
%   outputs. This is faster than requesting all the outputs.
%     Rather than providing a stride, you can alternatively specify a set of
%   columns. If that set of columns is just a single one, provide a third
%   argument, doStride=0, that tells Read to interpret the single number as a
%   column index rather than a stride.
%     After the third argument, one can pass the pair 'ridxs',value, where
%   value is a set of along-fault indices whose values to read. If not
%   passed, then all the along-fault cells are read.
%
% [ld c] = ReadConsts(datFn)
%   Read in the mat file. The first four ouptuts are parameters that are
%   passed to other functions in this file. ld is a structure containing
%   everything in the mat file. c is a struct of physical parameters.
%
% [mu mu_psi] = Calc_mu(c,is,psi,gamma)
%   Calculate mu and the partial derivative of mu w.r.t. psi.
%     gamma = chi - psi
%     is is an array of fault indices.
%     psi and gamma have length length(is).
%
% tau = Stress   (t,psi,chi,slip,p,c,[includeDamping])
%   New version of Stress to account for BEM stress calculations.
%  
% mr = MomentRate(x,v)
%   x is the along-fault cell-centered coordinates. v is slip rate formatted
%   as (along fault)x(time). Use the trapezoid rule to integrate.
  
  [varargout{1:nargout}] = feval(fn,varargin{:});
end
  
function [t psi chi dlt slip] =...
      Read(datFn,skip,doStride,varargin)
  % Input args
  ridxs = [];
  if(nargin > 3) ridxs = process_options(varargin,'ridxs',[]); end
  if(nargin < 2) skip = 1; end
  if(nargin < 3) doStride = length(skip) == 1; end
  % Consts
  [~,c] = ReadConsts(datFn,1);
  Ncell = c.Ncell;
  % Figure out which rows to read
  if(isempty(ridxs)) ridxs = 1:Ncell; end
  n = nargout;
  if(~isempty(ridxs))
    idxs = 1;
    for(i = 1:n-1)
      idxs = [idxs 1+(i-1)*Ncell+ridxs];
    end
  else
    idxs = 1:(1+(n-1)*Ncell);
  end
  % Mask against the number of rows available in case p,T data missing
  nr = SaveStreamData('NbrRows',datFn);
  idxs = intersect(1:nr,idxs);
  % Read
  A = SaveStreamData('Read',datFn,skip,idxs,doStride);
  % Extract
  t = A(1,:);
  if(nargout == 1) return; end
  lri = length(ridxs);
  idxs = 1:lri;
  psi  = A(1+      idxs,:);
  if(nargout == 2) return; end
  chi  = A(1+  lri+idxs,:);
  if(nargout == 3) return; end
  dlt = A(1+2*lri+idxs,:);
  slip = A(1+3*lri+idxs,:);
end
  
function [ld c] = ReadConsts(datFn,read_c)
  if(nargin < 2) read_c = 1; end
  if(read_c)
    ldfn = @load;
  else
    % Don't read c
    ldfn = @(fn)loadnot(fn,{'c'});
  end
  try
    [p n e] = fileparts(datFn);
    if(isempty(p)) p = '.'; end
    ld = ldfn([p filesep n '.mat']);
  catch
    ld = ldfn([datFn '.mat']);
  end
  c = [];
  if(~read_c) return; end
  c = ld.c;
  if(~isfield(c,'cbem'))
    if(isfield(ld.gP,'cbem'))
      c.cbem = ld.gP.cbem;
      % Zero out the multithreaded h structs
      c.cbem.h(2).ns = [];
      c.cbem.h(4).ns = [];
      c.cbem.thidx = 1;
    else
      c.cbem.use = 0;
    end
  end
  ld = rmfield(ld,'gP');
  if(c.cbem.use && isfield(c,'Gs'))
    c = rmfield(c,'Gs');
  end    
  ld.c = c;
end
  
function [tau mu] = Stress(t,psi,chi,slip,pf,c,incDamping,T)
  if(nargin < 8) T = []; end
  if(nargin < 7) incDamping = true; end
  n  = size(psi,2);
  v  = c.v0*exp(psi);
  mu = Calc_mu(c,1:c.Ncell,psi,chi-psi,T);
  s_normal = repmat(c.s_normal,1,n);
  if(~isempty(pf))
    s_eff = s_normal - pf;
  else
    s_eff = s_normal;
  end
  tau   = s_eff.*mu;
  if(incDamping)
    tau = tau + c.eta.*v;
  end
end
  
function tau = StressTransform(c,t,slip)
    s = [slip(:); c.v_creep*t];
    if (isfield(c,'delta_tau_fn') && ~isempty(c.delta_tau_fn))
      for n=1:length(t);
        des(:,n) = c.delta_tau_fn(c,t(n),0);
      end
    else des=0;
    end
    tau = c.tau0 + BEM_matvec(c,s,1) + des;
end

function [f f_psi f_gamma f_T] = Calc_mu(c,is,psi,gamma,T)
  deriv = nargout > 1;
  np = size(psi,2);
  a = repmat(c.a(is),1,np);
  b = repmat(c.b(is),1,np);
  if(length(c.mu_0) > 1)
    mu_0 = repmat(c.mu_0(is),1,np);
  else
    mu_0 = c.mu_0;
  end
  g = 0.5*exp(psi + (mu_0 + b.*gamma)./a);
  f = a.*asinh(g);
  if(deriv)
    gasinh_g = g./sqrt(g.^2 + 1);
    f_psi    = a.*gasinh_g;
    f_gamma  = gasinh_g.*b;
    ag_a = -g.*(mu_0 + b.*gamma)./a;
    f_a = asinhg + ag_a./den_asinh_g;
    a_T = a0/293;
    f_T = f_a.*a_T;
  end
end
  
function s = qload(fn,skip,doStride)
  if(nargin < 2) skip = 1; end
  if(nargin < 3) doStride = length(skip) == 1; end
  [s.t s.psi s.chi dlt s.slip] =...
      Read(fn,skip,doStride);
  s.fn = fn;
  s.gam = s.chi - s.psi;
  s.m = size(s.psi,1);
  s.n = length(s.t);
  [s.lo s.c] = ReadConsts(fn);
  s.lo = rmfield(s.lo,'c');
  s.es = repmat(s.c.s_normal,1,s.n);
  s.v = s.c.v0*exp(s.psi);
  s.x = CC(s.lo.x);
  if(length(s.x) > 1 && length(s.x) == size(s.psi,1))
    s.mr = MomentRate(s.x,s.v);
  elseif(length(s.x) == 1)
    s.mr = s.v;
  end
  [s.tau s.mu] = Stress(s.t,s.psi,s.chi,s.slip,[],s.c,1);
  if(length(skip) == 1)
    s.stride = skip;
    skip = 1:skip:skip*s.n;
  end
  s.datIdxs = skip;
  s.s2y = 1/(365.25*24*60*60);
  s.xkm = s.x*1e-3;

  [foo hn] = system('hostname');
  hn(hn < 'A' | hn > 'z') = [];
  s.hostname = hn;
end

function mr = MomentRate(x,v)
% x is the along-fault cell-centered coordinates. v is slip rate formatted as
% (along fault)x(time). Use the trapezoid rule to integrate.
  dx = diff(x(:));
  mr = sum((v(1:end-1,:) + v(2:end,:))/2 .* repmat(dx,1,size(v,2)))/...
       (x(end) - x(1));
end

function c = CC(v)
% Cell-centered from node-centered
  c = (v(1:end-1) + v(2:end))/2;
end

function ld = loadnot(fn,nlflds)
% Load a mat file, but exclude those variables listed in the cell array
% nlflds.
  flds = whos('-file',fn);
  flds = {flds.name};
  for(j = 1:length(nlflds))
    for(i = 1:length(flds))
      if(strcmpi(flds{i},nlflds{j}))
	flds(i) = [];
	break;
      end
    end
  end
  ld = load(fn,flds{:});
end

