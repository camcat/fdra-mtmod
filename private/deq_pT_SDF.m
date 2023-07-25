function stop = deq_pT_SDF(t,dt,y,msg,c,varargin)
% OutputFcn (see odeset). (t,y) are the new data.
  global gP;
  persistent etLastSave doDisp dispEvery saveEvery ...
      maxTimeBtwnSv profPts doDispText;
  tolfail = false;
  switch(msg)
   case 'init'
    % doDisp: Show plots
    % ssd:  SaveStreamData object
    % Drop last element of varargin, as that is extraneous stuff
    [doDisp dispEvery saveEvery maxTimeBtwnSv,...
      doDispText name] =...
	process_options(...
	    varargin(1:end-1),'disp',true,'dispEvery',50,...
	    'saveEvery',500,'maxTimeBtwnSv',1e4,...
	    'dispText',50,'name','');
    etLastSave = 0;
    gP.sdf.tf = t(end);
    return;
   case 'tolfail'
    tolfail = true;
   case 'done'
    return;
  end

  if(length(t) > 1)
    t = t(1);
    y = y(:,1);
  end
  
  slip  = y(1:c.Ncell);
  gamma = log(y(c.Ncell+1:2*c.Ncell));
  
  % If t gets too big, the step size can lose too many digits of
  % precision. >=v1.8: Now that we're using our own ode routine, I pass dt to
  % the ode function and so have an accurate dt for the diffusion equations
  % regardless of restarts. But still restart because the ode routine itself
  % maintains an hmin that is O(t).
  stop = 0;
  if(1e11*dt <= t || tolfail)
    gP.InitialStep = max(1e-10,dt);
    gP.y_init = [slip; exp(gamma); exp(gP.psi)];
    stop = 1;
  end
  if(tolfail) return; end
  
  doSaveNow = false;
  etLastSave = etLastSave + dt;
  if(etLastSave >= maxTimeBtwnSv || mod(gP.ctr,saveEvery) == 0 ||...
     t == gP.sdf.tf)
    etLastSave = 0;
    doSaveNow = true;
  end
  
  if(doSaveNow)
    chi = gamma + gP.psi;
    y = [gP.psi; chi; gamma; slip];
  end
  
  if(mod(gP.ctr,doDispText) == 0)
    fprintf(1,'% 7d %13.8e %13.6e',gP.ctr,gP.t_g+t,dt);
    % Time step stats
    prop = 0;
    fprintf(1,'%5.1f % 6.1f ',toc(gP.sdf.tic),prop);
    gP.sdf.tic = tic;
    % Accuracy
    fprintf(1,'%1.1e',gP.ode.mcon_psi);
    fprintf(1,'\n');
  end    

  if(doSaveNow)
    pp = []; pT = [];
    if(~isempty(profPts))
      pp = repmat(c.pinf(profPts)',c.sNy,1) + gP.p1(:,profPts);
      pT = repmat(c.Tinf(profPts)',c.sNy,1) + gP.T1(:,profPts);
    end
    data = [gP.t_g+t; y];
    gP.ssd.sv = SaveStreamData('Write',gP.ssd.sv,data);
    gP.nsv = gP.nsv + 1;
    
  end
    
  gP.t0 = t;
  gP.psi0 = gP.psi;
  gP.ctr = gP.ctr + 1;
  
  if(stop)
    % Issue an empty warning so that the while loop in fdra doesn't think
    % there was an integration failure.
    warning('');
  end
  
