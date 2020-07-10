function varargout = matrixfun(fun,varargin)
%MATRIXFUN Apply a function on a stack of matrices
% fun is a function handle
% varargin{i} contains the array Ai of size [nRows_i nColumns_i SZ]
% SZ = [sz1 sz2 ...] is [1] or a constant over the input matrices
% UniformOutput is true by default

% UniformOutput being set ?
uoArg = strcmp(varargin,'UniformOutput') ;
if any(uoArg) 
    uoArg = find(uoArg) ;
    UniformOutput = varargin(uoArg+1) ;
    varargin(uoArg+(0:1)) = [] ;
else
    UniformOutput = true ;
end

% Deal with sizes..
    SZ = cellfun(@size,varargin,'UniformOutput',false) ;
    NDIMS = cellfun(@numel,SZ) ;
% Take the biggest array as reference (not general...)
    [~,biggest] = max(NDIMS) ;
    SZ = SZ{biggest}(3:end) ;
    N = prod(SZ) ;


% Init Output
    varargout = cell(1,nargout) ; % Support multiple outputs
    out = varargout ;
% Apply to each matrix of the array
    in = cell(1,nargin-1) ;
    for nn = 1:N
        for ii = 1:nargin-1 ; in{ii} = varargin{ii}(:,:,nn) ; end
        [out{:}] = fun(in{:}) ;
        for oo = 1:nargout ; varargout{oo}{nn} = out{oo} ; end
    end

% Uniform output ?
    if UniformOutput
        for oo = 1:numel(varargout)
            varargout{oo} = reshape(cat(1,varargout{oo}{:}),[size(varargout{oo}{1}) SZ 1 1]) ;
        end
    else
        for oo = 1:numel(varargout)
            varargout{oo} = reshape(varargout{oo},[SZ 1 1]) ;
        end
    end

end

