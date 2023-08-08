function varargout = selectout(func,outputNo,varargin)
%function varargout = selectout(func,outputNo,varargin)
    varargout = cell(max(outputNo),1);
    [varargout{:}] = func(varargin{:});
    varargout = varargout(outputNo);
end