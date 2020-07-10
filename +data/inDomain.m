function [ival,iint] = inDomain(values,domains)
%IDOMAIN return indices of nd-values contained in nd-domains
% values: [nVal nCoord]
% domains [nInt 2 nCoord] (minValue_d maxValue_d)
%
% Alternative implementation to the following "brute-force" test:
%   dom = permute(domains,[2 1 3]) ; % [2 nInt nCoord]
%   INSIDE = all(dom(1,:,:)<=values,3) & all(dom(2,:,:)>=values,3) ;
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
iint = repelem(repmat(1:nInt,[1 nCoord]),nIndicesInInterval(:)) ; % lot of elements !

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

