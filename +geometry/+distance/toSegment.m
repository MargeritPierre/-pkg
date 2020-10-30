function d = toSegment(P,P1,P2)
%TOSEGMENT Distance from a list of points P to a (list of) segment(s) S = [P1->P2]
% inputs:
%   - P [nP nCoord]
%   - P1 [1 nCoord nS]
%   - P2 [1 nCoord nS]
% output: 
%   - distance d [nP 1 nS] 

    if isempty(P) ; d = [] ; return ; end

    [d,t] = pkg.geometry.distance.toLine(P,P1,P2) ;
    dp1 = pkg.geometry.distance.toPoint(P,P1) ;
    dp2 = pkg.geometry.distance.toPoint(P,P2) ;
    d(t<=0) = dp1(t<=0) ;
    d(t>=1) = dp2(t>=1) ;

end

