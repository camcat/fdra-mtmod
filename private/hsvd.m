function varargout = hsvd(fn,varargin)
% Compression of the BEM matrix to speed up MVPs.
%
% h = Build(A,err)
%   Construct an hsvd struct.
% h = RebuildThreaded(h,nthreads)
%   Rebuild the struct, separating it into pieces such that in the MVP y =
%   A*x, h is partitioned according to the partition for y.
% InitMex(h,id)
%   Initialize the mex structure. id = 0,1
% y = hsvd_mex(h.id,x)
% CleanupMex(id)
  [varargout{1:nargout}] = feval(fn,varargin{:});
end

function h = Build(A,err)
% We assume, as is true in FDRA, that the 1D BEM matrix is arranged spatially.
  h.sz = size(A);
  Amax = max(abs(A(:)));
  err = err*Amax/sqrt(prod(size(A)));
  n = NewLeaf([1 h.sz(1)],[1 h.sz(2)]);
  h.err = err;
  h.par = [1 h.sz(1)+1];
  % Decomposition
  szmin = 40;
  szmax = 700;
  h.ns = SplitRecurse(n,szmin,szmax);
  % Numerics
  for(i = 1:length(h.ns))
    n = h.ns(i);
    B = A(n.r(1):n.r(2),n.c(1):n.c(2));
    sz = size(B);
    [U S V] = svd(B,0);
    s = diag(S);
    if(false)
      % Find the number of singular values we need based on the Fro norm.    
      j = find(cumsum(s(end:-1:1).^2)/prod(size(B)) <= err^2,1,'last');
      j = length(s) - j;
    else
      j = maxerr(B,U,s,V',err);
    end
    if(sz(1)*j + sz(2)*j < sz(1)*sz(2))
      h.ns(i).U = U(:,1:j)*diag(s(1:j));
      h.ns(i).Vt = V(:,1:j)';
    else
      h.ns(i).B = B;
    end
  end
  h.ns = {h.ns};
end

function j = maxerr(A,U,s,Vt,err)
% Return the minimum j such that
%   max_{pixel} |U(:,1:j)*diag(s(1:j))*V(:,1:j)' - A| <= err.
  B = zeros(size(A));
  j = 0;
  while(j < length(s))
    if(max(abs(A(:)-B(:))) <= err) return; end
    j = j + 1;
    B = B + s(j)*U(:,j)*Vt(j,:);
  end
end

function h1 = RebuildThreaded(h,A,nthreads)
  % Partition the rows
  np = round(h.sz(1)/nthreads);
  par = (0:nthreads-1)*np + 1;
  par(end+1) = h.sz(1) + 1;
  % Break up blocks that cross a partition boundary
  h1 = h;
  h1.par = par;
  h1.ns = cell(1,nthreads);
  for(i = 1:length(h.ns))
    for(j = 1:length(h.ns{i}))
      n = h.ns{i}(j);
      p1 = find(n.r(1) >= par,1,'last');
      p2 = find(n.r(2) >= par,1,'last');
      if(p1 == p2)
	n.tid = p1;
	% Same partition, so do nothing
	h1.ns{n.tid} = [h1.ns{n.tid}; n];
      else
	for(k = 1:p2-p1+1)
	  n1 = Partition(n,[par(p1+k-1) par(p1+k)-1],A);
	  n1.tid = p1 + k - 1;
	  h1.ns{n1.tid} = [h1.ns{n1.tid}; n1];
	end
      end
    end
  end
end

function n = Partition(n,par,A)
  % New rows
  r1 = n.r(1);
  if(n.r(1) < par(1)) n.r(1) = par(1); end
  if(n.r(2) > par(2)) n.r(2) = par(2); end 
  rs = n.r - r1 + 1;
  % Extract part of the block
  if(isempty(n.B))
    k = size(n.Vt,1);
    sz = [n.r(2)-n.r(1)+1 n.c(2)-n.c(1)+1];
    if(sz(1)*k + sz(2)*k < sz(1)*sz(2))
      n.U = n.U(rs(1):rs(2),:);
    else
      % Block is now small enough that it's better to use the original
      % block
      n.B = A(n.r(1):n.r(2),n.c(1):n.c(2));
      n.U = []; n.Vt = [];
      if(n.r(2) < n.r(1))keyboard;end
    end
  else
    n.B = n.B(rs(1):rs(2),:);
  end
end

function InitMex(h,id)
% Initialize the mex file for performing MVPs.
  hsvd_mex(h,'init',id);
end

function CleanupMex(id)
  hsvd_mex(0,'cleanup',id);
end

function ns = SplitRecurse(n,szmin,szmax)
  % Base case
  if(n.r(2) - n.r(1) + 1 < 2*szmin || n.c(2) - n.c(1) + 1 < 2*szmin)
    ns = n;
    return;
  end
  % Recurse
  rh = Half(n.r);
  ch = Half(n.c);
  ns = [NewLeaf([n.r(1)   rh],[n.c(1)   ch]) % (1,1)
	NewLeaf([rh+1 n.r(2)],[ch+1 n.c(2)]) % (2,2)
	NewLeaf([n.r(1)   rh],[ch+1 n.c(2)]) % (1,2)
	NewLeaf([rh+1 n.r(2)],[n.c(1) ch])]; % (2,1)
  ns = [SplitRecurse(ns(1),szmin,szmax)
	SplitRecurse(ns(2),szmin,szmax)
	SplitIfTooBig(ns(3),szmax)
	SplitIfTooBig(ns(4),szmax)];
end

function ns = SplitIfTooBig(n,szmax)
  % Base case
  if(n.r(2) - n.r(1) + 1 <= szmax && n.c(2) - n.c(1) + 1 <= szmax)
    ns = n;
    return;
  end
  % Recurse
  rh = Half(n.r);
  ch = Half(n.c);
  ns = [NewLeaf([n.r(1)   rh],[n.c(1)   ch]) % (1,1)
	NewLeaf([rh+1 n.r(2)],[ch+1 n.c(2)]) % (2,2)
	NewLeaf([n.r(1)   rh],[ch+1 n.c(2)]) % (1,2)
	NewLeaf([rh+1 n.r(2)],[n.c(1) ch])]; % (2,1)
  ns = [SplitIfTooBig(ns(1),szmax)
	SplitIfTooBig(ns(2),szmax)
	SplitIfTooBig(ns(3),szmax)
	SplitIfTooBig(ns(4),szmax)];
end

function n = NewLeaf(r,c);
  n = struct('r',r,'c',c,'U',[],'Vt',[],'B',[],'tid',1);
end

function h = Half(r)
  h = round((r(1) + r(2))/2);
end

function M = Show(h,A)
  subplot(221); imagesc(log10(abs(A))); axis image; axis off;
  title('log_{10}|A|');
  
  M = [];
  for(i = 1:length(h.ns))
    for(j = 1:length(h.ns{i}))
      n = h.ns{i}(j);
      M(n.r(1):n.r(2),n.c(1):n.c(2)) = rand;
    end
  end
  subplot(222); imagesc(M); axis image; axis off;
  title('Blocks');
  
  M = zeros(size(A));
  Mth = M;
  for(i = 1:length(h.ns))
    for(j = 1:length(h.ns{i}))
      n = h.ns{i}(j);
      if(isempty(n.B))
	B = n.U*n.Vt;
      else
	B = n.B;
      end
      M(n.r(1):n.r(2),n.c(1):n.c(2)) = B;
      Mth(n.r(1):n.r(2),n.c(1):n.c(2)) = n.tid;
    end
  end
  subplot(223); imagesc(log10(abs(A-M)./max(abs(A(:)))));
  axis image; axis off; colorbar;
  ne = nnz(h);
  title(sprintf('Relative error with %3.1fx compression',...
		prod(size(A))/ne));
  subplot(224); imagesc(Mth); axis image; axis off;
  title('Threads');
end

function ne = nnz(h)
  ne = 0;
  for(i = 1:length(h.ns))
    for(j = 1:length(h.ns{i}))
      n = h.ns{i}(j);
      ne = ne + prod(size(n.U)) + prod(size(n.Vt)) + prod(size(n.B));
    end
  end
end

function y = mvp(h,x)
% Reference Matlab implementation
  y = zeros(h.sz(1),size(x,2));
  for(i = 1:length(h.ns))
    for(j = 1:length(h.ns{i}))
      n = h.ns{i}(j);
      if(isempty(n.B))
	upd = n.U*(n.Vt*x(n.c(1):n.c(2),:));
      else
	upd = n.B*x(n.c(1):n.c(2),:);
      end
      y(n.r(1):n.r(2),:) = y(n.r(1):n.r(2),:) + upd;
    end
  end
end

function [yt yh] = TestMex(A,h)
  h.id = 0;
  InitMex(h,h.id);
  x = randn(size(A,2),1);
  tic;
  for(i = 1:100) yt = A*x; end
  etA = toc/100;
  fprintf(1,'A: %1.3e\n',etA);
  tic;
  for(i = 1:100) yh = hsvd_mex(h.id,x); end
  etm = toc/100;
  fprintf(1,'h: %1.3e\n',etm);
  fprintf('speedup = %3.1f\n',etA/etm);
  CleanupMex(h);
  ym = mvp(h,x);
  semilogy(abs((yt-yh)./yt)); axis tight;
  title(sprintf('Ax vs h(x). compress=%3.1f speedup=%3.1f',...
		prod(size(A))/nnz(h),etA/etm));
end
