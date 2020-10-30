function d = toPoint(P,PT)
%TOPOINT Distance from a list of points P to a (list of) point(s) PT
% inputs:
%   - P [nP nCoord]
%   - PT [1 nCoord nPT]
% output: 
%   - distance d [nP 1 nL] 

    if isempty(P) ; d = [] ; return ; end

    d = sqrt(sum((P-PT).^2,2)) ;

end

