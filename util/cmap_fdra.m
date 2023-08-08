function cm = cmap_fdra(v);
%function cm = cmap_fdra([vels]);
%vels contains: [Vlo Vcr Vss Vco Vfast] = locked velocity (blue), creeping velocity (green), slow slip events (yellow), seismic velocity (red), fast seismic velocity (dark red).
%default 10.^[-16, -9, -6, -1, 0].

if nargin<1 v=10.^[-16, -9, -6, -1, 0]'; end

Vlo=log10(v(1));
Vpl=log10(v(2));
Vss=log10(v(3));
Vco=log10(v(4));
Vf=log10(v(5));

Ntot=1e3;

N1=int16((Vpl-Vlo)/(Vco-Vlo)*Ntot);
N2=int16((Vss-Vpl)/(Vco-Vlo)*Ntot);
N3=int16((Vco-Vss)/(Vco-Vlo)*Ntot);
N4=int16((Vf-Vco)/(Vco-Vlo)*Ntot);

cm1=viridis(int16(N1*6/5));
cm1=cm1(1:N1,:);
cm2=viridis(int16(N2*6));
cm2=cm2(end-N2:end,:);

r=[1 0 0]; %red
dr=[0.2 0 0]; %dark red

cm3=interp1(double([0 N3]), [cm2(end,:); r], double(1:1:N3));
cm4=interp1(double([0 N4]), [r; dr], double(1:1:N4));

cm=[cm1; cm2; cm3; cm4];

