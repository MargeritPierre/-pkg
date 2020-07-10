function [ival,iint] = inDomain(values,domains)
%IDOMAIN return indices of nd-values contained in nd-domains
% values: [nVal nCoord]
% domains [nInt 2 nCoord] (minValue_d maxValue_d)

% Infos
nVal = size(values,1) ;
nInt  = size(domains,1) ;
nCoord = min(size(values,2),size(domains,3)) ;

% Format
values = values(:,1:nCoord) ;
domains = domains(:,1:2,1:nCoord) ;

% Put values and interval bounds together
data = [values ; reshape(domains,[2*nInt nCoord])] ; % [nVal+2*nInt nCoord]
% Keep indices corresponding to interval bounds
ibounds = reshape((1:2*nInt)+nVal,[nInt 2]) ; % [nInt 2]

% Sort data
[~,isort] = sort(data,1) ;
% Retrieve indices corresponding to values
isVal = isort<=nVal ;

% Number of values before each index
np = cumsum(isVal(:)) ;
% Sorted index of each interval bound
ii = NaN(nInt,2,nCoord) ;
for cc = 1:nCoord ; [~,ii(:,:,cc)] = ismember(ibounds,isort(:,cc)) ; end
ii = ii + reshape((2*nInt+nVal)*(0:nCoord-1),[1 1 nCoord]) ;
% Number of values before the interval bound
npi = np(ii) ; % [nInt 2 nCoord]
npi = permute(npi,[1 3 2]) ; % [nInt nCoord 2]

% Long vector of interval indices
nIndicesInInterval = diff(npi,1,3) ;
iint = repelem(repmat(1:nInt,[1 nCoord]),nIndicesInInterval(:)) ;

% Long vector of value indices
starts = npi(:,:,1)+1 ; ends = npi(:,:,2) ;
lengths = ends(:)'-starts(:)'+1 ; cumLengths = cumsum(lengths) ;
vv = (0:cumLengths(end)-1) + repelem(starts(:)' - [0 cumLengths(1:end-1)],lengths) ;

% Retrieve the corresponding value indices
ival = reshape(isort(isVal),[nVal nCoord]) ;
ival = ival(vv) ;

if nCoord>1 % multiple coordinates
    % Accumulation sparse array
    A = sparse(ival(:),iint(:),1,nVal,nInt) ;
    % Find values being in intervals for all coordinates
    [ival,iint,nValidCoord] = find(A) ;
    valid = nValidCoord==nCoord ;
    ival = ival(valid) ; iint = iint(valid) ;
end

end

