function [ival,iint] = inDomain(values,domains,method)
%INDOMAIN return indices of nd-values contained in nd-domains
% values: [nVal nCoord]
% domains [nInt 2 nCoord] (minValue_d maxValue_d)
%
% Alternative implementation to the following "brute-force" test:
%   dom = permute(domains,[2 1 3]) ; % [2 nInt nCoord]
%   pt = reshape(values,[nVal 1 nCoord]) ; % [nVal 1 nCoord]
%   INSIDE = all(dom(1,:,:)<=pt,3) & all(dom(2,:,:)>=pt,3) ;
%	[ip,ie] = find(INSIDE) ;
%
% It is memory-efficient and the speedup depends on the sparsity of INSIDE
%   high sparsity: up to two order of magnitude faster!
%   low sparsity: hem.. not so efficient!
%

% Infos
nVal = size(values,1) ;
nInt  = size(domains,1) ;
nCoord = min(size(values,2),size(domains,3)) ;

% Format
if size(values,2)>nCoord ; values = values(:,1:nCoord) ; end
if size(domains,3)>nCoord ; domains = domains(:,1:2,1:nCoord) ; end

if nargin<3 ; method = 'sort' ; end
switch method
    case 'sort'
        coordSorting ;
    case 'tree'
        jointTree ;
    otherwise
        bruteForce ;
end


% END OF MAIN FUNCTION
% ===================================================


%% NESTED FUNCTIONS 


% BRUTE FORCE ALGORITHM: TEST ALL vs. ALL
function bruteForce
    dom = permute(domains,[2 1 3]) ; % [2 nInt nCoord]
    pt = reshape(values,[nVal 1 nCoord]) ; % [nVal 1 nCoord]
    INSIDE = all(dom(1,:,:)<=pt,3) & all(dom(2,:,:)>=pt,3) ;
    [ival,iint] = find(INSIDE) ;
end


% FIND POINTS IN DOMAIN BY FIRST SORTING (MEMORY EFFICIENT)
% critical point: array accumulation (sparse)
function coordSorting

    % Put values and interval bounds together
    data = [values ; reshape(domains,[2*nInt nCoord])] ; % [nVal+2*nInt nCoord]
    nData = nVal+2*nInt ;
    % Keep indices corresponding to interval bounds
    ibounds = reshape(nVal+1:nData,[nInt 2]) ; % [nInt 2]

    % Sort data
    [~,isort] = sort(data,1) ; % [nData nCoord]
    % Retrieve indices corresponding to values
    isVal = isort<=nVal ; % [nData nCoord]

    % Number of values before each index
    np = cumsum(isVal(:)) ; % [nData*nCoord 1]
    % Sorted index of each interval bound
    ii = NaN(nInt,2,nCoord) ;
    for cc = 1:nCoord ; [~,ii(:,:,cc)] = ismember(ibounds,isort(:,cc)) ; end
    ii = ii + reshape((2*nInt+nVal)*(0:nCoord-1),[1 1 nCoord]) ; % [nInt 2 nCoord]
    % Number of values before the interval bound
    npi = np(ii) ; % [nInt 2 nCoord]
    npi = permute(npi,[1 3 2]) ; % [nInt nCoord 2]

    % Long vector of interval indices
    nIndicesInInterval = diff(npi,1,3) ; % [nInt nCoord]
    iint = repelem(repmat(1:nInt,[1 nCoord]),nIndicesInInterval(:)) ; % a lot of elements !

    % Long vector of value indices
    starts = npi(:,:,1)+1 ; ends = npi(:,:,2) ; % [nInt nCoord]
    lengths = ends(:)'-starts(:)'+1 ; cumLengths = cumsum(lengths) ; % [nInt*nCoord]
    vv = (0:cumLengths(end)-1) + repelem(starts(:)' - [0 cumLengths(1:end-1)],lengths) ; % lot of elements ! (critical assignment)

    % Retrieve the corresponding value indices
    ival = isort(isVal) ; % [nVal*nCoord]
    ival = ival(vv) ; % lot of elements ! (critical assignment)

    if nCoord>1 % multiple coordinates
        % Accumulation sparse array
        A = sparse(ival(:),iint(:),1,nVal,nInt) ; % (critical)
        % Find values being in intervals for all coordinates
        [ival,iint] = find(A==nCoord) ;
    end

end


% USE TEST DECIMATION BY JOINT R-TREE SPLITTING
function jointTree
    
    % Reshape data
    dom = permute(domains,[2 1 3]) ; % [2 nInt nCoord]
    pt = reshape(values,[nVal 1 nCoord]) ; % [nVal 1 nCoord]

    % Put values and interval bounds together
    data = [values ; reshape(domains,[2*nInt nCoord])] ; % [nVal+2*nInt nCoord]
    nData = nVal+2*nInt ;
    % Keep indices corresponding to values and intervals
    %iInt = [NaN(1,nVal) 1:nInt 1:nInt] ;
    %iVal = [1:nVal NaN(1,2*nInt)] ;
    
    % Joint group splitting (median)
    % groups will contain values and one ore two interval bounds
    N = 3*nCoord ; % final number of groups = 2^N
    groups = {1:nData} ; % start with the full group
    for nn = 1:N
        coord = mod(nn-1,nCoord)+1 ; % coordinate for the median
        splitGroups = {} ; % init new groups
        for g = 1:numel(groups)
            ind = groups{g} ; % select the indices in the group
            med = median(data(ind,coord),1) ; % split with median value
            higher = (data(ind,coord)>=med)' ; % coordinates higher than median
            isUpperBnd = ind>nVal+nInt ; % is an interval upper bound
            isLowerBnd = ind>nVal & ~isUpperBnd ; % is an interval lower bound
            upperPart = higher | (isLowerBnd & ~higher & ~ismember(ind+nInt,ind(isUpperBnd))) ;
            lowerPart = (~higher) | (isUpperBnd & higher & ~ismember(ind-nInt,ind(isLowerBnd))) ;
            splitGroups = [splitGroups ind(lowerPart) ind(upperPart)] ; % fill groups
        end
        groups = splitGroups ;
    end

    % Process groups
    nGroups = numel(groups) ;
    for gg = 1:nGroups
        ii = groups{gg} ;
        isValue = ii<=nVal ;
        ival = ii(isValue) ;
        iint = mod(ii(~isValue)-nVal-1,nInt)+1 ;
        iint = unique(iint) ; %numel(iint)
        inside = all(dom(1,iint,:)<=pt(ival,:,:),3) & all(dom(2,iint,:)>=pt(ival,:,:),3) ;
        [vvv,iii] = find(inside) ;
        groups{gg} = [reshape(ival(vvv),[],1) reshape(iint(iii),[],1)] ;
    end
    groups = cat(1,groups{:}) ;
    ival = groups(:,1) ;
    iint = groups(:,2) ;
    
    
    
    
end

end
% END OF CODE
% ===================================================


function tests
%% TESTS ===================================================
clearvars
nC = 2 ;
nP = 1000000 ;
nB = 10000 ; Bsz = 0.01 ;

P = rand(nP,nC) ;
B = rand(nB,1,nC) + Bsz*[-1 1].*rand(nB,1,nC) ;

if nC>1 && nC<4
    clf
    axis equal
    patch('vertices',P,'faces',(1:nP)','marker','.','MarkerEdgeColor','r')
    ind = permute(B,[1 3 2]) ;
    ind = [ind(:,:,1) ; [ind(:,1,2) ind(:,2,1)] ; ind(:,:,2) ; [ind(:,1,1) ind(:,2,2)]] ;
    patch('vertices',ind,'faces',(1:nB)'+(0:3)*nB,'FaceColor','none','EdgeColor','b')
end

clc
profile on
%tic ; [ival1,iint1] = pkg.data.inDomain(P,B,'brute') ; toc
%tic ; [ival1,iint1] = pkg.data.inDomain(P,B,'sort') ; toc
tic ; [ival2,iint2] = pkg.data.inDomain(P,B,'tree') ; toc
profile off

%sparse(ival1,iint1,true,nP,nB)-sparse(ival2,iint2,true,nP,nB)



%% PLOT GROUPS
clf ;
axis equal
for gg = 1:numel(groups)
    ii = groups{gg} ;
    isValue = ii<=nVal ;
    ival = ii(isValue) ;
    med = plot(values(ival,1),values(ival,2),'.') ;
    iint = mod(ii(~isValue)-nVal-1,nInt)+1 ;
    iint = unique(iint) ;
    ind = permute(domains(iint,:,:),[1 3 2]) ;
    ind = [ind(:,:,1) ; [ind(:,1,2) ind(:,2,1)] ; ind(:,:,2) ; [ind(:,1,1) ind(:,2,2)]] ;
    patch('vertices',ind,'faces',(1:numel(iint))'+(0:3)*numel(iint),'FaceColor','none','EdgeColor',med.Color)
end

%% SPARSITY INDICATORS
minOvlp = reshape(sum(diff(B,1,2),1),[],nC)./range(reshape(B,[],nC),1) ;
sparsity = sum(prod(diff(B,1,2),3))/(prod(range(reshape(B,[],nC),1),2)*nB) ;





%% END OF TESTS
end





