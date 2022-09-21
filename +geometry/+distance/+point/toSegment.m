function [d,t] = toSegment(P,P1,P2)
%TOSEGMENT Distance from a list of points P to a (list of) segment(s) S = [P1->P2]
% inputs:
%   - P [nP nCoord]
%   - P1 [(nP or 1) nCoord nS]
%   - P2 [(nP or 1) nCoord nS]
% output: 
%   - distance d [nP 1 nS] 
%   - parameter p [nP 1 nS] 

    if isempty(P) ; d = [] ; return ; end

    [d,t] = pkg.geometry.distance.point.toLine(P,P1,P2) ;
    dp1 = pkg.geometry.distance.point.toPoint(P,P1) ;
    dp2 = pkg.geometry.distance.point.toPoint(P,P2) ;
    d(t<=0) = dp1(t<=0) ;
    d(t>=1) = dp2(t>=1) ;
    t = min(max(t,0),1) ;

end

