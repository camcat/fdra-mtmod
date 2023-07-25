function fdra(sFn,varargin)
% FDRA1D 1.33: One-Dimensional Fault Dynamics with a Radiation Damping
%              Approximation
%
% This is a minimal version modified from Andrew Bradley's original code.
%
% fdra(scriptFn,{'name',val}*)
%
% The parameters are as follows:
%
% Sizes:
%   N_cell: x direction
%
% Software:
%   saveFn: save files prefix
%   allow_overwrite: Allow saveFn to be overwritten?
%   conType = 'loose', 'moderate', 'tight'
%   save_every: write frequency in times steps
%   maxTimeBtwnSv: approx upper bound, in seconds, between saved time steps
%   dispText: display text? if >0, display every dispText steps.
%   use_nthreads: number of OpenMP CPU threads
%
% Initial values:
%   psi_init, chi_init, p_init, T_init, optionally slip_init
%   ts, tend
%     p_init and T_init are total values, including T_inf.
%  
% Physics:
%  General scalars:
%   eta, V0, G, poisson, dip
%  General vectors:
%   Dc, s_normal
%
%  Rate-and-state friction:
%   evolution = 1 (aging), 2 (slip)
%   mu_0, a, b
%
%  Elasticity (Half-space BEM):
%   x: (N_cell+1)-vector of cell boundaries
%   use_compressed_bem: use a faster BEM computation?
%   Vpl_creep: downdip loading velocity
%
% fdra produces three output files. Let saveFn='fn'. Then the output files
% are
%    fn.mat: Contains all the variables defined in script, as well as several
%      more. The useful ones include
%        - c: A structure of physical constants
%          - c.map{k} is the set of fault node indices to which sd*.op(k)
%            applies. These sets of nodes are determined automatically based
%            on fault physical parameter values provided in the user script.
%        - script: A string containing a copy of the user script file used
%          for these results.
%        - options: ode23 options struct.
%    fn.dat: Fault variables saved every save_every time steps. See
%      Util('Read') to learn how to extract these.
%
% Authors:
%   A.M. Bradley (ambrad@cs.stanford.edu)
%   P. Segall    (segall@stanford.edu)
%   C. Cattania  (camcat@mit.edu)
%   with help from S. Schmitt, J. Chen, L. Bruhat
%   CDFM, Dept of Geophysics, Stanford University
  
  %Shorthand for finding cell centers
  CC = @(x) (x(2:end)+x(1:end-1))/2;

  % If the mex files are not present, build them
%  BuildMex();
  
  fprintf(1,'... Processing user input.\n');
    
  saveFn = sFn.saveFn;
  if ~isfield(sFn,'allow_overwrite') sFn.allow_overwrite = false; end
  if (~sFn.allow_overwrite && exist([saveFn '.mat'],'file'))
    error(sprintf('File %s exists, my friend!\n', saveFn));
  end
    
  % Pack constants
  c.Ncell     = length(sFn.d_c);
  c.eta       = sFn.eta;
  c.v0        = sFn.V0;
  c.evolution = sFn.evolution;
  c.mu_0      = sFn.mu_0(:);
  c.d_c       = sFn.d_c(:);
  c.s_normal  = sFn.s_normal(:);
  c.bes       = c.s_normal;
  c.a         = sFn.a(:);
  c.b         = sFn.b(:);
  c.G 	      = sFn.G;
  c.poisson   = sFn.poisson;
  c.v_creep   = sFn.Vpl_creep;

  x = sFn.x;
  ts = sFn.ts;
  tend = sFn.tend;
  psi_init = sFn.psi_init(:);
  chi_init = sFn.chi_init(:);

  %by default load both sides and fullspace, and use plane strain conditions on a dipping fault
  if isfield(sFn,'dip') c.dip = sFn.dip; else dip=10; end
  if isfield(sFn,'fullspace') c.fullspace = sFn.fullspace; else c.fullspace = true; end
  if isfield(sFn,'antiplane') c.antiplane = sFn.antiplane; else c.antiplane = false; end
  if isfield(sFn,'use_nthreads') use_nthreads = sFn.use_nthreads; else use_nthreads=0; end
  if isfield(sFn,'use_compressed_bem') use_compressed_bem = sFn.use_compressed_bem; 
     else use_compressed_bem = 0; end
  if isfield(sFn,'load_bothsides') c.load_bothsides = sFn.load_bothsides;
     else c.load_bothsides = true; end

  % disp: Plot results during integration. dispEvery: Display every dispEvery steps.
  % saveEvery, maxTimeBtwnSv: Save every saveEvery steps or when
  %   maxTimeBtwnSv simulation seconds have elapsed. maxTimeBtwnSv is only an
  %   approximate upper bound on the time between steps.
  if isfield(sFn,'maxTimeBtwnSv') maxTimeBtwnSv = sFn.maxTimeBtwnSv; else maxTimeBtwnSv = 6e5; end
  if isfield(sFn,'dispText') dispText = sFn.dispText; else dispText = 1; end
  if isfield(sFn,'save_every') save_every = sFn.save_every; else save_every = 100; end

  clear global gP;
  global gP;
  
  % Convergence parameters
  % TODO reinstate these
%  fprintf(1,'Tolerance: %s\n',conType);
%  gP.conType = lower(conType);
%  [c.PsiTol c.PsiMaxIts c.pTTol c.RelTol c.hmTol] = Tols(conType);
%  c.fail_tol = 1e-7;
  
%   if (BEM_warn(xb)) return; end

  c.Index = 1:c.Ncell;

  xb = x(c.Index);
  xb(end+1) = x(c.Index(end)+1);    
  c.x=xb;
  c.Gs = GetTractions2d(xb, CC(xb), c.dip, c.G, c.poisson,1,c.antiplane, c.fullspace, c.load_bothsides);
  c.cbem = CompressedBEM(c,use_compressed_bem,use_nthreads);
  
  gP.t_g = ts;
  gP.InitialStep = 2.5e-5;
  gP.nsv = 0;
  gP.bes = c.bes;
  gP.ctr = 0;
  gP.psi = psi_init;
  gP.psi0 = gP.psi;
  
  odeFn = @deq_pT;

  % Flag to let deq_pT know it's in a fail state. Reset when yset is requested.
  gP.ode.fail = 0;
  gP.ode.mcon_psi = 0;
  
  % Initial values
  slip_init = zeros(c.Ncell,1);
  % Get tau0, the initial stress that gives psi
  c.tau0 = Util('Stress',ts,psi_init,chi_init,slip_init,[],c,1);
  gamma_init = chi_init - psi_init;
  gP.y_init = [slip_init(:); exp(gamma_init(:)); exp(psi_init(:))];
  y_nonneg = ones(size(gP.y_init));
  abstol = 1e-40*ones(size(gP.y_init));
  
  % [t psi chi dlt slip]
  gP.ssd.sv = SaveStreamData('Init',[saveFn '.dat']);
  
  % Clean up memory
  if (use_compressed_bem)
    % Clear Gs from memory
    c = srmfield(c,'Gs');
  end

  o = ones(c.Ncell,1);
  %options = odeset(...
  %    'RelTol',c.RelTol,'AbsTol',abstol,...
  %    'NonNegative',find(y_nonneg));
  save([saveFn '.mat'],'-v7.3');
  
  fprintf(1,'Saving to %s\n',saveFn);
  fprintf(1,['   Step           Time   Step length  npslv   CPU     mx ',...
	     'err-psi  err-pT\n']);

  gP.sdf.tic = tic;
  intfail = 0;
  while (true)
%    options.InitialStep = gP.InitialStep;
    gP.t0 = 0;   

%TODO set options for matlab ode23.

    options.OutputFcn = @(t,y,flag) deq_pT_SDF_matlab (t,y,flag,c,...
          'disp',true,'dispEvery',save_every,...
          'dispText',dispText,'saveEvery',save_every,...
          'maxTimeBtwnSv',maxTimeBtwnSv,...
          'name',saveFn,varargin);

    ops.OutputFcn=options.OutputFcn;
    ode23(@(t,y) odeFn(t,y,c),[gP.t0 tend - gP.t_g],gP.y_init,ops);

    gP.t_g = gP.t_g + gP.t0;
    if (gP.t_g >= tend) break; end
    % Integration failure?
    if (~isempty(strfind(lastwarn,'Failure at')))
      intfail = intfail + 1;
    else
      intfail = 0;
      if (c.print_debug >= 2)
	fprintf(1,'Restarting for time step tol.\n');
      end
    end
    if (intfail >= 1)
      if (gP.ctr <= 2)
	msg = [...
  'This can happen at the start of a new simulation if a phsyical\n',...
  'parameter has a value that is greatly different than what we have up\n',...
  'to now considered reasonable.'];
      else
	msg = [...
  'All went well until now. Did something very unstable just happen?'];
      end
      fprintf(1,sprintf(...
	  'Quitting because of integration failure.\n  %s\n',msg));
      break;
    end
  end
end

function cbem = CompressedBEM(c,use_compressed_bem,use_nthreads)
  cbem.use = use_compressed_bem;
  if (use_compressed_bem)
    error('Error: use_compressed_bem not implemented.');
  end
%    fprintf(1,'... Compressing GF arrays.\n');
%    bemstr = 'Compressed ';
%    err = 1e-9; % down from 1e-12 after change in hsvd.m
%    cbem.h = BuildCBEM(c.Gs,err,[0 1],use_nthreads);
%    cbem.strt = [0 2]; % start-1 of the Gs h structs
%    % thread index
%    if (use_nthreads <= 1)
%      cbem.thidx = 1; 
%    elseif (~c.no_pT)
%      % Should be 2 in principle. But I found that multithreaded hsvd does
%      % not interact well with multithreaded deq_solve_*.c, and the latter is
%      % much more time consuming. Hence make hsvd always singlethreaded.
%      cbem.thidx = 1;
%    else
%      % p,T are not being evolved, so we might as well multithread hsvd.
%      cbem.thidx = 2;
%    end
%    if (cbem.thidx == 1)
%      cbem.h(2).ns = [];
%      cbem.h(4).ns = [];
%    end
%  else
%    bemstr = '';
%  end
%  fprintf(1,'%sBEM.\n',bemstr);
end

% ------------------------------------------------------------------------------
% BEM matrix

function [Gs Gn] = GetTractions2d (xb, xo, dip, mu, nu, want_vbc, antiplane, fullspace, load_bothsides)
  if (nargin < 6) want_vbc = true; end

  if (antiplane)
    if (~fullspace & abs(dip-90)>10*eps) error('Screw dislocation on dipping fault not implemented.'); end
    t = repmat(xo(:) , 1, numel(xb));
    z = repmat(xb(:)', numel(xo), 1);
    
    if (fullspace) Gs = mu/(2*pi) * [1./(t-z)];  %TODO check this
       else Gs = mu/(2*pi) * [1./(t-z) - 1./(t+z)];
    end

  else
    dip = dip*pi/180;
    cd  = cos(  dip);
    sd  = sin(  dip);
    c2d = cos(2*dip);
    c4d = cos(4*dip);
  
    t = repmat(xo(:) , 1, numel(xb));
    z = repmat(xb(:)', numel(xo), 1);
 
    if (fullspace)
      Gs = (mu/(1-nu))/(2*pi) * [1./(t-z)];
    else 
      Gs = (-2.*z.*(t + z).*mu.*...
         (t.^4 - 3.*t.^3.*z + 7.*t.^2.*z.^2 -...
          3.*t.*z.^3 + z.^4 + (t.^4 - 5.*t.^3.*z + 4.*t.^2.*z.^2 -...
                               5.*t.*z.^3 + z.^4).*...
          c2d + t.^2.*z.^2.*c4d).*sd.^2)./...
        (pi.*(t - z).*(-1 + nu).*(t.^2 + z.^2 - 2.*t.*z.*c2d).^3);  
    end
  end

  if load_bothsides
     Gload=Gs(:,1)-Gs(:,end);
  else
    Gload=-Gs(:,end); 
  end

  Gs = [-Gs(:,1:end-1) + Gs(:,2:end), Gload];
  if (~want_vbc) Gs = Gs(:,1:end-1); end
  
  if (nargout > 1)
    if (antiplane)
      Gn = 0*Gs;
    else
    Gn = (4.*z.^2.*mu.*cd.*(-3.*t.^2.*z + z.^3 + 2.*t.^3.*c2d).*sd.^3)./...
         (pi.*(-1 + nu).*(t.^2 + z.^2 - 2.*t.*z.*c2d).^3);
    if load_bothsides
       Gload=Gn(:,1)-Gn(:,end);
    else
       Gload=-Gn(:,end); 
    end
    Gn = [-Gn(:,1:end-1) + Gn(:,2:end), Gload];
    if (~want_vbc) Gn = Gn(:,1:end-1); end
  end
end
end
