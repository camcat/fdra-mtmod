%A couple of settings and shortcuts for convenience

year2sec = 365.25*24*60*60;
sec2year = 1/year2sec;
cc = @(x) (x(2:end)+x(1:end-1))/2;    %cell centers
Linf = @(G, nu, dc, a, b, sigma) G/(1-nu)*dc/pi./sigma.*(b./(b-a).^2);
Lc = @(G, nu, dc, b, sigma) 1.377*G/(1-nu)*dc./(b.*sigma);

findpos = @(x, pos)  selectout(@(x) min(x), 2, abs(x-pos));
lineyk = @(y) plot(get(gca,'xlim'), [y y], 'k','HandleVisibility','off'); 
lineyd = @(y) plot(get(gca,'xlim'), [y y], '--k','HandleVisibility','off'); 
linexk = @(x) plot([x x], get(gca,'ylim'),'k','HandleVisibility','off'); 
linexd = @(x) plot([x x], get(gca,'ylim'),'--k','HandleVisibility','off'); 

idx = @(q,x) findpos(q.xkm, x);  %give position in km
idt = @(q,t) findpos(q.t*sec2year, t);	 %give time in yrs

qload = @(fn) Util('qload', fn);

