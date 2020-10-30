function [indices,starts,ends] = linspace(starts,ends,lengths)
%LINSPACE Linspace operator accepting vector inputs
% examples <TODO>: 
% - [indices,starts,ends] = pkg.data.linspace([0 5 3],[3 2 7],[1 -1 2])
%       -indices: [0 1 2 3 5 4 3 2 3 5 7]
%       -starts: [1 5 9]
%       -ends: [4 8 11]
% - [indices,starts,ends] = pkg.data.linspace([0.1 3],[0.5 0],[0.1 -1.5])
%       -indices: [.1 .2 .3 .4 .5 3 1.5 0]
%       -starts: [1 6]
%       -ends: [5 8]

% Input checks
if nargin<3 ; error('not enough arguments') ;end
if nargin<3 ; lengths = ones(size(starts))*100 ; end

% Format inputs
O = (starts(:)+ends(:)+lengths(:))'*0 ;
starts = starts(:)' + O ;
ends = ends(:)' + O ;
lengths = lengths(:)' + O ;

% Individual lengths
lengths = floor(lengths) ;
lengths(lengths<0) = 0 ; % negative lengths->empty vectors
cumLengths = cumsum(lengths) ; % cummuative lengths

% Steps
steps = (ends-starts)./(lengths-1) ;

% Indices
indices = 0:cumLengths(end)-1 ; % original indices
indices = indices - repelem([0 cumLengths(1:end-1)],lengths) ; % shift to zero
indices = indices .* repelem(steps,lengths) ; % apply steps
indices = indices + repelem(starts,lengths) ; % shift starting indices

% starts & ends output
starts = [0 cumLengths(1:end-1)] + 1 ;
ends = cumLengths ;

end

