function ls=plotasp(q);
%function ls=plotasp(q)

shorthands

pts=find(sign(abs(diff(q.c.a-q.c.b)))>0);
dx=mean(diff(q.xkm)); %Will assume it's uniform;

hold on
for n=1:length(pts)
  ls{n}=linexk(q.xkm(pts(n))+0.5*dx);
  set(ls{n},'LineWidth',1);
end
