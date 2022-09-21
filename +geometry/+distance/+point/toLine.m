function [d,t] = toLine(P,P1,P2)
%TOLINE Distance from a list of points P to a (list of) infinite line(s) L = [P1->P2]
% inputs:
%   - P [nP nCoord]
%   - P1 [(nP or 1) nCoord nL]
%   - P2 [(nP or 1) nCoord nL]
% outputs: 
%   - distance d [nP 1 nL] 
%   - normalized parameter t [nP 1 nL] 

    if isempty(P) ; d = [] ; return ; end

    u = (P2-P1)./sqrt(sum((P2-P1).^2,2)) ;
    t = sum((P-P1).*u,2) ;
    d = sqrt(abs(sum((P-P1).^2,2) - t.^2)) ;
    if nargout>1 ; t = t./sqrt(sum((P2-P1).^2,2)) ; end

end

