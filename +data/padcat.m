function [DATA,id] = padcat(dim,varargin)
%PADCAT concatenate an inhomogeneous set of data 
% size are first homogeneized
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
data = varargin ;

% Case where data is not a cell
if ~iscell(data) 
    id = 1 ; 
    return ;
end

% Keep data chunks indices
data = data(:) ;
id = 1:numel(data) ;

% Sizes
sz = cellfun(@size,data,'UniformOutput',false) ;
nd = cellfun(@ndims,data) ;

% Homogenize size vector
finalDim = max(dim,max(nd)) ; 
% Mask for valid size elements
    mask = zeros(finalDim,numel(data)) ;
    mask(sub2ind(size(mask),nd(:)',1:numel(data))) = 1 ;
    mask = cumsum(mask,1,'reverse') ;
    mask = logical(mask) ;
% Size matrix
    SZ = ones(finalDim,numel(data)) ;
    SZ(mask) = cat(2,sz{:}) ;
    SZ = SZ' ;

% Keep indices of data chunks
id = repelem(id(:),SZ(:,dim),1) ;

% <TODO> Can we regularly concatenate the data ?
%[uSZ,ia,ic] = unique(SZ,'rows','stable') ;

% Initialize the DATA
szData = max(SZ,[],1) ; 
szData(dim) = sum(SZ(:,dim)) ;
DATA = repmat(filler,szData) ;

% nd-Mask
MASK = zeros(szData) ;
lastInd = repelem(SZ,SZ(:,dim),1) ;
lastInd(:,dim) = 1:size(lastInd,1) ;
lastInd = num2cell(lastInd,1) ;
MASK(sub2ind(szData,lastInd{:})) = 1 ;
for dd = [1:dim-1 dim+1:ndims(MASK)]
    MASK = cumsum(MASK,dd,'reverse') ;
end
MASK = logical(MASK) ;

% Reshaping
data = cellfun(@reshape,data,repmat({[]},numel(data),1),repmat({1},numel(data),1),'UniformOutput',false) ;
data = cat(1,data{:}) ;
DATA(MASK) = data ;

end

