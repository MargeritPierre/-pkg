function varargout = matchsize(dims,varargin)
%MATCHSIZE match the sizes of an inhomogeneous set of data 
% if 'filler' is in varargin, then the next argument is the (scalar) filler

% Filler ?
isfiller = strcmp(varargin,'filler') ; 
if any(isfiller)
    fillIdx = find(isfiller) ;
    filler = varargin{fillIdx+1} ;
    varargin(fillIdx+[0 1]) = [] ;
else
    filler = NaN ;
end

% Case where on cell only have been given
if numel(varargin)==1 && iscell(varargin) 
    varargin = varargin{1} ;
end

% Check outputs
if nargout>1 && nargout~=numel(varargin)
    error('Wrong output format.') ;
end

% PROCESS
if ~isempty(dims)
    
    % Sizes
    sz = cellfun(@size,varargin,'UniformOutput',false) ;
    nd = cellfun(@ndims,varargin) ;

    % Homogenize size vector
    finalDim = max(max(dims),max(nd)) ; 
    % Mask for valid size elements
        mask = zeros(finalDim,numel(varargin)) ;
        mask(sub2ind(size(mask),nd(:)',1:numel(varargin))) = 1 ;
        mask = cumsum(mask,1,'reverse') ;
        mask = logical(mask) ;
    % Size matrix
        SZ = ones(finalDim,numel(varargin)) ;
        SZ(mask) = cat(2,sz{:}) ;
        SZ = SZ' ;

    % Homogenized target sizes
    hSZ = SZ ;
    hSZ(:,dims) = repmat(max(SZ(:,dims),[],1),[size(SZ,1) 1]) ;

    % Homogenize
    for ii = 1:numel(varargin)
        varargin{ii} = padarray(varargin{ii},hSZ(ii,:)-SZ(ii,:),filler,'post') ;
    end

end

% Output
if nargout<=1
    varargout = {varargin} ;
else
    varargout = varargin ;
end

end

