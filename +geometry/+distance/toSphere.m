function d = toSphere(P,CENTER,RADIUS)
%TOSPHERE SIGNED Distance from a list of points P to a (list of) Sphere(s)
% inputs:
%   - P [nP nCoord]
%   - CENTER [1 nCoord nC] center of the sphere(s)
%   - RADIUS [1 nCoord nC] radius of the sphere(s)
% output: 
%   - d [nP 1 nC] signed distance (d<0 inside, d>0 outside)

    if isempty(P) ; d = [] ; return ; end

    d = pkg.geometry.distance.toPoint(P,CENTER)-RADIUS ;

end

