function [indices,starts,ends] = colon(starts,ends,steps)
%COLON Colon operator accepting vector inputs
% examples: 
% - [indices,starts,ends] = pkg.data.colon([0 5 3],[3 2 7],[1 -1 2])
%       -indices: [0 1 2 3 5 4 3 2 3 5 7]
%       -starts: [1 5 9]
%       -ends: [4 8 11]
% - [indices,starts,ends] = pkg.data.colon([0.1 3],[0.5 0],[0.1 -1.5])
%       -indices: [.1 .2 .3 .4 .5 3 1.5 0]
%       -starts: [1 6]
%       -ends: [5 8]

% Input checks
if nargin<2 ; error('not enough arguments') ;end
if nargin<3 ; steps = ones(size(starts)) ; end

% Empty case
if isempty(starts)
    indices = [] ; starts = [] ; ends = [] ;
    return ;
end

% Individual lengths
lengths = floor((ends-starts)./steps)+1 ;
lengths(lengths<0) = 0 ; % negative lengths->empty vectors
cumLengths = cumsum(lengths) ; % cummuative lengths

% Indices
indices = 0:cumLengths(end)-1 ; % original indices
indices = indices - repelem([0 cumLengths(1:end-1)],lengths) ; % shift to zero
if nargin>2 ; indices = indices .* repelem(steps,lengths) ; end % apply steps
indices = indices + repelem(starts,lengths) ; % shift starting indices

% starts & ends output
starts = [0 cumLengths(1:end-1)] + 1 ;
ends = cumLengths ;

end

