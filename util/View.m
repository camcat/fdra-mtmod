function [fignos ] = myView(s , ops, flds, idx, idt);
% function [ ] = myView( s, ops, flds, idx, idt)
% makes colorplots of various quantities
% if multifig==1, makes separate figures; otherwise, a single figure with
% subplots.
%tau_c: see Paul's SCEC proposal.

y2s = 365.25*24*60*60;
s2y = 1/y2s;

if (nargin<2)
    multifig=1;
    truespace=1;
    truetime=0;
    rescale_v=0;    
    tau_c=0;
    force_highres=0;    %plot all time steps
    originalsteps=0;
    make_subplot=0;
    logtimeT0=0;
    timeunit='y';
    switchxy=0;
else
    if ~isfield (ops,'switchxy') switchxy=0;          
    else switchxy=ops.switchxy; end
    if ~isfield (ops,'logtime') logtimeT0=0;          
    else logtimeT0=ops.logtime; end
    if ~isfield (ops,'timetofailure') timetofailure=1; %for backwards compatilibity.
    else timetofailure=ops.timetofailure; end
   %NB originalsteps get overridden by trueXX
    if ~isfield (ops,'originalsteps') originalsteps=0;          
    else originalsteps=ops.originalsteps; end
    if ~isfield (ops,'multifig') multifig=1; 
    else multifig=ops.multifig; end
    if ~isfield (ops,'truespace') truespace=1;
    else truespace=ops.truespace; end
    if ~isfield (ops,'timeunit') timeu='y';
    else timeu=ops.timeunit; end
    if ~isfield (ops,'truetime') truetime=0;
    else truetime=ops.truetime; end
    if ~isfield (ops,'v_over_v0') rescale_v=0;
    else rescale_v=ops.v_over_v0; end
    if ~isfield (ops,'tau_c') tau_c=0; 
    else tau_c=ops.tau_c; end
    if ~isfield (ops,'sel_ts') 
    else
        s=Util_camcat('SelectTimes',s,ops.sel_ts); 
    end
    if ~isfield (ops,'force_highres') force_highres=0; 
    else force_highres=ops.force_highres; end
    if ~isfield (ops,'subplot') make_subplot=0;          
    else make_subplot=ops.subplot; end
end

if (std(diff(s.x))/mean(diff(s.x)) > 1e-8) var_mesh=1;
else var_mesh=0;
end

if (nargin<3 || isempty(flds)) 
    flds=struct('v',{});
end

if (nargin<4 || isempty(idx)) idx=1:1:length(s.xkm); end
if (nargin<5 || isempty(idt)) 
	if logtimeT0 
          if timetofailure idt=find(s.t<logtimeT0);
          else idt=find(s.t>logtimeT0); end 
	else idt=1:1:length(s.t); end
end

figure
cmap=colormap(cmap_fdra(10.^[-12.5, -7.5, -3, -1, 1]));
close

    function makecplot (x, y0, z0)
        if (multi_file) idx0=fidx;
        else idx0=ones(size(y0));
        end        
        N=length(unique(idx0));
        %idx contains indices of elements for each subplot.
        for n=1:N
          y=y0(idx0==n);
          z=z0(idx0==n,:);
   	  if switchxy
	      xt=x;
	      x=y;
	      y=xt;
	      z=z';
	  end
          if(~make_subplot) subplot(N,1,N-n+1); end
          if (force_highres) 
	      disp('Forcing high resolution');
 	      if (truetime | (truespace & var_mesh) | logtimeT0)  pcolor(x,y,z);
	      else imagesc(x,y,z); set(gca,'YDir','normal'); 
	      end
          else  mypcolor(x,y,z);
          end
          shading flat; 
          colormap(cmap);
          cb=colorbar('Location','NorthOutside');
	  if (logtimeT0 & timetofailure)
	     if switchxy set(gca,'XDir','reverse');
	     else set(gca,'YDir','reverse');
	     end
	  end
          if (n==N)
              if truespace xlab='Distance (km)';
              else xlab='Distance (grid points)'; end
          end
          if truetime ylab = ['Time (', timeu, ')'];
	  elseif logtimeT0 ylab = 'Log_{10} Time';
	  else ylab = 'Time steps'; end

	  if switchxy xlabel(ylab); ylabel(xlab);
	  else ylabel(ylab); xlabel(xlab);
	  end
        end
    end

if (~multifig) subplot (2,2,1); end

multi_file=0;

if (truetime) 
 switch (timeu)
 case 's'
	ts=s.t(idt); 
 case 'y'
	ts=s.t(idt)*s2y;
 otherwise 
 	error('Invalid time unit: must be s or y.');
 end
elseif (logtimeT0)
 if timetofailure
   idt(s.t(idt)>logtimeT0)=[];
   ts=log10(logtimeT0-s.t(idt));
 else 
   idt(s.t(idt)<logtimeT0)=[];
   ts=log10(-logtimeT0+s.t(idt));
 end
else
    if (originalsteps)
        ts=s.datIdxs(idt);
        fidx=s.fileIdx(idt);
        if (length(unique(s.fileIdx(idt)))>1)
            multi_file=1;
            if (make_subplot)
                error('Multiple files, option originalsteps=1 and multi_subplot=1not compatible.')    
            else
                warning('Multiple files and option originalsteps=1: will make subplots')                
            end
        end
    else
        ts=idt;
    end
end
if (truespace) xs=s.xkm(idx); 
else
    if (originalsteps)
        try
            xs=s.ridxs(idx);
        catch
            xs=idx;   %variable does not exist if loaded 'qload' insead of 'qloadr';
        end
    else
        xs=idx;
    end
end
  
allfld=fieldnames(flds);
for n=1:length(allfld)  

    if (make_subplot==0)
        figure
        gg=gcf;
        fignos(n)=gg.Number;
    else
        subplot(make_subplot(1), make_subplot(2), make_subplot(3:end))
    end
    
    %velocity:    
    if (strcmp(allfld{n},'v'))
        if rescale_v
         makecplot(xs, ts, log10(abs(s.v(idx,idt))'/s.lo.Vpl_creep));
         caxis([-10 0]/s.lo.Vpl_creep);
        else
         makecplot(xs, ts, log10(abs(s.v(idx,idt))'));
         caxis([-14 0]);
         cb.Title.String='log10(V, m/s)';
        end        
    %stress:
    elseif (strcmp(allfld{n},'tau'))
        if (tau_c)
            %steady state value for plate creep velocity (with time dep. normal
            %stress)
            tss=repmat(f_ss(s),1,length(s.t)).*s.es;
            makecplot(xs, ts, log10(s.tau(idx,idt)'./tss'));
            title('\tau/\tau_{ss.creep}')
        end    
        %tau0s=repmat(s.c.tau0,1,size(s.tau,2));
        %makecplot(xs, ts, s.tau'-tau0s');
        %title('\tau-\tau_0')
        makecplot(xs, ts, s.tau(idx,idt)');
        title('\tau')        
   %strength:
   elseif (strcmp(allfld{n},'strength'))
        tau=s.es.*s.mu;
        makecplot(xs, ts, tau(idx,idt)');
        title('\tau_f = \mu\sigma_{eff}')        
    else
        try
           sallfldn=s.(allfld{n});        
        catch   %own field, not found in s.
           sallfldn=flds.(allfld{n});  
        end        
        makecplot(xs, ts, sallfldn(idx,idt)');
       title(allfld{n})       
        
    end
    
     try
        caxis(flds.(allfld{n}));
     end

end

end

function [argout] = mypcolor (varargin)
%Will use pcolor if the matrix is small, pcolor_small otherwise.
%arguments must be in the form mypcolor(X,Y,Z, [...]).

%number of time steps:
nts=length(varargin{2});
nts_max=500;    %number of time steps for plotting.

if (nts)>nts_max
    pcolor_small(ceil(nts/nts_max),1,varargin{:});
else
    pcolor(varargin{:})
end

end
