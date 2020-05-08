function [pp,dist] = closestPoint(PQ,P,method)
%CLOSESTPOINT Return the indices of the closest points P to PQ
% Chooses the fastest method between direct vectorization, loop over P or
% loop over PQ
%
% INPUTS
% P = [nP nCoord]
% PQ = [nPQ nCoord]
% method = -1 (auto choice), 
%          0 (test all), 
%          1 (vectorization), 
%          2 (loop over PQ), 
%          3 (loop over P) 
%
% OUTPUTS
% pp = [nPQ 1] indices such that dist(PQ(ii),P(pp(ii))) = min(dist(PQ(ii),P)
% dist = [nPQ 1] corresponding distance

% Method choice parameters (see below)
maxArraySize = 1e8 ;

% Sizes
nP = size(P,1) ;
nPQ = size(PQ,1) ;
nCoord = numel(P)/nP ;
arraySize = nP*nPQ*nCoord ;

% Method choice
if arraySize<=maxArraySize ; choice = 1 ; 
elseif 1 || nPQ<=nP ; choice = 2 ; % Seems to be faster than method 3 most of the cases
else ; choice = 3 ;
end
% Apply if needed
if nargin<3 ; method = -1 ; end
if method == -1 ; method = choice ; end

% Timing
test = method==0 ;
if test ; disp('') ; disp(['closestPoint - timing (choice: ' num2str(choice) ')']) ; end


% Direct Vectorization
if method == 1 || (test && arraySize<=maxArraySize)
    if test ; tic ; end
    
    PQ = reshape(PQ,[nPQ 1 nCoord]) ; % [nPQ 1 nCoord]
    P = reshape(P,[1 nP nCoord]) ; % [1 nP nCoord]
    [dist,pp] = min(sum((PQ-P).^2,3),[],2) ; % [nPQ 1]
    
    if test ; disp(['  1|vectorized: ' num2str(toc*1000,3) ' ms']) ; 
    else ; return ; end
end

% Initialize for loops
dist = NaN(nPQ,1) ;
pp = ones(nPQ,1) ;
PQ = reshape(PQ,[nPQ nCoord]) ; % [nPQ nCoord]
P = reshape(P,[nP nCoord]) ; % [nP nCoord]

% Loop over PQ
if method == 2 || test
    if test ; tic ; end
    
    for ii = 1:nPQ
        [dist(ii),pp(ii)] = min(sum((PQ(ii,:)-P).^2,2)) ;
    end
    
    if test ; disp(['  2|loop.PQ: ' num2str(toc*1000,3) ' ms']) ; 
    else ; return ; end
end

% Loop over P
if method == 1 || test
    if test ; tic ; end
    
    for ii = 1:nP
        newDist = sum((PQ-P(ii,:)).^2,2) ;
        jj = find(newDist<dist) ;
        pp(jj) = ii ;
        dist(jj) = newDist(jj) ; 
    end
    
    if test ; disp(['  3|loop.P: ' num2str(toc*1000,3) ' ms']) ; 
    else ; return ; end
end

