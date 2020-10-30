function d = toCircle(P,CENTER,RADIUS)
%TOCIRCLE SIGNED Distance from a list of points P to a (list of) Circle(s)
% inputs:
%   - P [nP nCoord]
%   - CENTER [1 nCoord nC] center of the circle(s)
%   - RADIUS [1 nCoord nC] radius of the circle(s)
% output: 
%   - d [nP 1 nC] signed distance (d<0 inside, d>0 outside)

    if isempty(P) ; d = [] ; return ; end

    d = pkg.geometry.distance.toPoint(P,CENTER)-RADIUS ;

end

