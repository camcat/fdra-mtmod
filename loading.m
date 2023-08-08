function out = loading(c,t,deriv,mode)
%This is a function defining the stressing rates. Its arguments are:
%c: structure containing grid points and other physical properties (produced by the input script)
%t: time
%deriv: boolean value. If deriv=0, returns stress; otherwise, returns stressing rate.
%
%the output is a scalar or an array of dimension size(c.x), with the stress or stressing rates at each grid point. Here it's simply set to constant rate.
%mode=1 (shear), 2 (normal).

if exist('mode')~=1 mode=1; end

switch (mode)

case 1
  taudot=c.taudot;
case 2
  taudot=c.sigmadot;
otherwise
  error('loading.m: mode must be 1 or 2');
end

if (deriv)
  out = taudot;
else 
  out = taudot*t;
end

