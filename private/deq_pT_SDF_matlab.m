function status = deq_pT_SDF_matlab(t,y,flag,varargin)
%function status = deq_pT_SDF_matlab(t,y,flag,varargin)
%save as deq_pT_SDF, but compatible with matlab's ode23.

if isempty(flag)
   deq_pT_SDF(t,0, y,'', varargin{:});
else
   deq_pT_SDF(t,0,y,flag, varargin{:});
end

status = 0;

