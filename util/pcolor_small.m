function [ varargout ] = pcolor_small( fact, dim, varargin)
%function [ varargout ] = pcolor_small(fact, dim, varargin)
% like pcolor, but reduces matrix sixe by a factor of 'fact' before
% plotting.
% if dim==1 or 2, will reduce along first (second)  dimension; if dim==12,
% will reduce along both (by an amount close to sqrt(fact) in each
% dimension).
% All standard pcolor variables can be appended.

%check if there are at least 3 input variables, and the first 3 are
%numeric
if (length(varargin)>=3 && sum(cellfun( @(x) isnumeric(x), varargin(1:3)))==3)
    M=varargin{3};
    Mpos1= false;
else
    M=varargin{1};
    Mpos1= true;
end

switch(dim)
    case 1 %reduce along first dimension
        M=M(1:fact:end,:);
        %reduce coordinate vector
        if (~Mpos1) varargin{2}=varargin{2}(1:fact:end); end
            
    case 2 %reduce along first dimension
        M=M(:,1:fact:end);
        if (~Mpos1) varargin{1}=varargin{1}(1:fact:end); end
    case 12 %reduce along both dimensions; it's a bit approximate.        
        f1=ceil(sqrt(fact));
        f2=floor(sqrt(fact));
        
        if (length(M(:,1)) > length(M(1,:)))     %first dimension larger
            M=M(1:f1:end, 1:f2:end);
            if (~Mpos1) varargin{1}=varargin{1}(1:f2:end); end
            if (~Mpos1) varargin{2}=varargin{2}(1:f1:end); end
        else
            M=M(1:f2:end, 1:f1:end);
            if (~Mpos1) varargin{1}=varargin{1}(1:f1:end); end
            if (~Mpos1) varargin{2}=varargin{2}(1:f2:end); end
        end                   
    otherwise
        error('invalid value for dim: should be 1, 2 or 12(to reduce along both directions)')        
        return
end

if (Mpos1) varargin{1} = M;
else varargin{3} = M;
end

if (nargout) varargout{1:nargout} = pcolor(varargin{:});
else  
	pcolor(varargin{:});
end

end

