function [eqks_start, eqks_end] = find_eqks(s, v_thres, mindt);
%function [eq0, eq1, nucl] = find_eqks(s, [v_thres=0.1], [mindt=1]);
%Detects periods in which somewhere along the fault the slip velocity
%exceeds v_thresh. Returns start, end indices of blocks for which this is
%true.

if nargin<2 v_thres=0.1; end
if nargin<3 mindt=1; end

eqks_blocks=max(s.v)>v_thres;
eqks_start=find(diff(eqks_blocks)==1)+1;
eqks_end=find(diff(eqks_blocks)==-1);

tooshort=1;
while ~isempty(tooshort)
  tooshort=find(diff(s.t(eqks_start))<mindt);
  eqks_start(tooshort+1)=[];
  eqks_end(tooshort)=[];
end

try 
 if eqks_start(1)>eqks_end(1) eqks_end(1)=[]; end
 if eqks_start(end)>eqks_end(end) eqks_start(end)=[]; end
end
