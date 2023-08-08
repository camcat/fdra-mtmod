function ls=plotrups(q, vthres)
%function ls=plotrups(q, vthres)
%Adds horizontal lines indicating seismic ruptures to color plot.

if nargin<2 vthres=0.1; end
shorthands

%Find earthquakes:
[es ee]=find_eqks(q);

hold on
for n=1:length(es)
  rup = find(max(q.v(:,es(n):ee(n))'>vthres));
  ls{n}=plot(q.xkm(rup([1 end])), q.t(ee(n))*[1 1]*sec2year,'-r', 'LineWidth',2);
end


