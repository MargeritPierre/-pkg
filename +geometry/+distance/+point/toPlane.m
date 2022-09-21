function d = toPlane(P,O,N)
%TOPLANE SIGNED Distance from a list of points P to a (list of) infinite plane(s)
% inputs:
%   - P [nP nCoord] % points
%   - O [1 nCoord nL] % plane origin(s)
%   - N [1 nCoord nL] % plane normal(s)
% outputs: 
%   - distance d [nP 1 nL] 
%   - signed distance sd [nP 1 nL]

    if isempty(P) ; d = [] ; return ; end

    N = N./sqrt(sum(N.^2,2)) ; % normalize normals
    d = sum((P-O).*N,2) ; % signed distance sd = (P-O).N

end

