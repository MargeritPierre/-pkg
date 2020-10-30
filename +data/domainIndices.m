function [ind,dom] = domainIndices(domains)
%INDOMAIN return nd-indices (integer valued) contained in nd-domains
% domains [nInt 2 nCoord] (minValue_d maxValue_d)

% Infos
nInt  = size(domains,1) ;

% Domain bounds
xmin = ceil(domains(:,1,:)) ; % [nInt 1 nCoord]
xmax = floor(domains(:,end,:)) ; % [nInt 1 nCoord]

% Number of indices in each domain, along each dimension
nIndInDomainDim = xmax-xmin+1 ; % [nInt 1 2]
% Number of indices in each domain
nIndInDom = prod(nIndInDomainDim,3) ; % [nInt 1 1]
% Total number of indices
nIndTotal = sum(nIndInDom,1) ; % [1 1 1]

% Keep domain indices
dom = repelem(1:nInt,nIndInDom)' ; % [nIndTotal 1]

% Concatenate all indices
nIndInPrevDomains = [0 cumsum(nIndInDom(1:end-1)')] ;
pp = (0:nIndTotal-1) - repelem(nIndInPrevDomains,nIndInDom) ; % linear index

% Subindices
cumSZ = [ones(nInt,1) cumprod(nIndInDomainDim(:,:),2)] ; % [nInt nCoord+1]
ind = mod(pp(:),repelem(cumSZ(:,2:end),nIndInDom,1)) ; % [nIndTotal nCoord]
ind = floor(ind./repelem(cumSZ(:,1:end-1),nIndInDom,1)) ; % [nIndTotal nCoord]
ind = ind + xmin(dom,:) ; % [nIndTotal nCoord]



end
% END OF CODE
% ===================================================


function tests
%% TESTS ===================================================
clearvars
nCoord = 2 ;
nDom = 100000 ; maxDomSz = 5 ; 1000^(1/nCoord) ;

DOM = round(maxDomSz*rand(nDom,1,nCoord) + 0.5*maxDomSz*[-1 1].*rand(nDom,1,nCoord)) ;


profile on
tic ; [ind,dom] = pkg.data.domainIndices(DOM) ; toc
profile off

numel(dom)

%% END OF TESTS
end





