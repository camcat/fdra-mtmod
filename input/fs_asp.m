function s = fs_asp(model);
% Input script to fdra. 

% These are some useful shortcuts
shorthands

% Path to output files:
homef='output/';
fname=['asp' num2str(model)];

% Evolution law: aging (1) slip (2)
s.evolution = 1;

% Constitutive parameters
s.mu_0      = 0.60;             % nominal friction coefficient
Vs	    = 3.464e3;		% shear wave speed
rho	    = 2670;		% density
s.G 	    = Vs^2*rho;		% shear modulus
s.Vpl_creep = 0;		% plate rate
s.V0	    = 1e-6;		% reference slip rate
s.poisson = 0.25;                % Poisson ratio 
s.eta     = 0.5*s.G/Vs;     % radiation damping parameter
s.dip       = nan;		% fault dip (not needed in fullspace)
Vinit	 = 1e-9; 	      % initial velocity
s.taudot = 3e-4; 		%loading (shear stressing rate)

s.antiplane=0;
s.fullspace=1;


if mod(model,2)==0
  Wf=40;
  Ar=10;
  Ax=20;
 %lc/dx, i.e. number of grid points in cohesive zone. 4-8 is reasonable.
  lcodx=6;
else
  Wf=120;		% Total fault length (km)
  Ar=[7 18 1];	% Asperity radii (km)
  Ax=[25 70 110];	% Asperity center position (km)
  lcodx=4;
end

% Frictional parameters 
b0=0.015;
avw=0.009;
sigma=5e7;          
d_c=3e-2; %Large value to keep computation time reasonable.


if model<2
  % Vel. strengthening patches
  avs=0.021;    
  lc=Lc(s.G, s.poisson, d_c, b0, sigma);
  dx=lc/lcodx;
  dxk=dx*1e-3; %km

  % Cell edges:
  x_flt=[0 dxk:dxk:Wf];
else
  sigma_low=5e6;
  
  % Start with low res mesh
  lc=Lc(s.G, s.poisson, d_c, b0, sigma_low);
  lca=Lc(s.G, s.poisson, d_c, b0, sigma);
  %dx=min(lc/lcodx, Wf*1e3/8);
  dx=lc/lcodx;
  dxk=dx*1e-3; %km
  dxa=lca/lcodx;
  dxa=dx * 2.^floor(log2(dxa/dx));
  dxak=dxa*1e-3;
  % Cell edges:
  x_flt=[0 dxak:dxak:Wf];
  
  % Remove points from stable regions until desired resolution is reached:
  while 2*max(diff(x_flt))<min(dxk,2)
      isasp=0*x_flt;
      for n=1:length(Ar)
         ii = find(x_flt>(Ax(n)-Ar(n)-dxak) & x_flt<(Ax(n)+Ar(n)+dxak));
         if isempty(ii)
             [~,ii]=min(abs(x_flt-Ax(n)))
         end
         i0 = max(ii(1)-1,1);
         i1 = min(ii(end)+1, length(x_flt));
         isasp(i0:i1)=1;
      end
      notasp=find(~isasp);
      x_flt(notasp(2:2:end))=[];
  end
end

dim=[1 length(x_flt)-1];
s.b=b0*ones(dim);
s.d_c=d_c*ones(dim);

if model<2 s.a=avs*ones(dim);
else s.a=avw*ones(dim);
end
if model<2 s.s_normal=sigma*ones(dim);
else s.s_normal=sigma_low*ones(dim);
end

for n=1:length(Ar)
  ii=x_flt>(Ax(n)-Ar(n)) & x_flt<(Ax(n)+Ar(n));
  if model<2 s.a(ii)=avw;
  else s.s_normal(ii)=sigma;
  end
end

x_flt=x_flt*1e3;   %km2m/:90
s.x = x_flt;
Nflt = length(s.x) - 1;

one = ones(1,Nflt);
s.psi_init = one*log(Vinit/s.V0);
s.chi_init = 2*one;

s.ts   = 0;
%End times gives few cycles with 3MPa stress drop. Difference between
%models (0,1) and (2,3) was found empirically.
s.tend = (1+double(model>1))*3e6/s.taudot;

%Setup up loading: constant shear stressing rate
s.delta_tau_fn= @(c,t,deriv) loading(c,t,deriv);


%--------------------------------------------------------------% 
%--------------------------------------------------------------%

s.use_compressed_bem = 1; % Compressed BEM (small error, large speed up)
s.save_every = 5;       % How often to save data in solver time steps
s.dispText = 100;         % How often to show output
s.conType = 'loose';      % Convergence tolerance: loose, moderate, tight

s.saveFn = [homef fname];
s.allow_overwrite = false; % Allow the file saveFn to be overwritten?

