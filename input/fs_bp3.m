function s = fs_bp3;
% Input script to fdra. 
% Parameters set as described here: https://strike.scec.org/cvws/seas/download/SEAS_BP3.pdf

% Path to output files:
homef='output/';
fname='test';

% Set default values:
dz=25;	%cell size (m)
psi = 90;  %dip angle (deg)

% Evolution law: aging (1) slip (2)
s.evolution = 1;
% Convenience
year2sec = 365.25*24*60*60;
sec2year = 1/year2sec; 
CC = @(x) (x(2:end)+x(1:end-1))/2;    %cell centers

% Constitutive parameters
s.mu_0      = 0.60;             % nominal friction coefficient
Vs	    = 3.464e3;		% shear wave speed
rho	    = 2670;		% density
s.G 	    = Vs^2*rho;		% shear modulus
s.Vpl_creep = 1e-9;		% plate rate
s.V0	    = 1e-6;		% reference slip rate
s.poisson = 0.25;                % Poisson ratio 
s.eta     = 0.5*s.G/Vs;     % radiation damping parameter
s.dip       = psi;
Vinit	 = 1e-9; 	      % initial velocity

s.antiplane=0;
s.fullspace=0;

% Frictional parameters 
b0=0.015;
amax=0.025;
a0=0.010;
sigma=5e7;          
d_c=8e-3;

Wf=40;	%km
H=15;	%km
h=3;	%km
dzk=dz*1e-3; %km

% Cell edges:
x_flt=[0 dzk:dzk:Wf];

dim=[1 length(x_flt)-1];
s.b=b0*ones(dim);
s.a=a0+(amax-a0)*(CC(x_flt)-H)/h;
s.a(CC(x_flt)<H)=a0;
s.a(CC(x_flt)>H+h)=amax;

s.d_c=d_c*ones(dim);
s.s_normal=sigma*ones(dim);
   
x_flt=x_flt*1e3;   %km2m
  
s.s_normal = s.s_normal;

s.x = x_flt;
Nflt = length(s.x) - 1;

one = ones(Nflt,1);
     
% Initial conditions
tau0 = sigma*amax * asinh( Vinit/(2*s.V0)*exp((s.mu_0 + b0*log(s.V0/Vinit))/amax) ) + s.eta*Vinit;
gamma_init = (s.a./s.b).*log(2*s.V0/Vinit*sinh((tau0-s.eta*Vinit)./(s.a*sigma))) - s.mu_0./s.b;

one = ones(1,Nflt);
s.psi_init = one*log(Vinit/s.V0);
s.chi_init = gamma_init+s.psi_init;

s.ts   = 0;
s.tend = 300*year2sec;  


%--------------------------------------------------------------% 
%--------------------------------------------------------------%

s.use_compressed_bem = 0; % Compressed BEM (small error, large speed up)
s.save_every = 5;       % How often to save data in solver time steps
s.dispText = 100;         % How often to show output
s.conType = 'loose';      % Convergence tolerance: loose, moderate, tight

s.saveFn = [homef fname];
s.allow_overwrite = false; % Allow the file saveFn to be overwritten?

