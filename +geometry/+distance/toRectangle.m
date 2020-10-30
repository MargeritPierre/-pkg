function d = toRectangle(P,CENTER,SIDES)
%TORECTANGLE SIGNED Distance from a list of points P to a Rectangle
% inputs:
%   - P [nP nCoord]
%   - CENTER [1 nCoord 1] center of the rectangle
%   - SIDES [1 nCoord 1] side dimensions of the rectangle
% output: 
%   - d [nP 1 1] signed distance (d<0 inside, d>0 outside)

    if isempty(P) ; d = [] ; return ; end

    PTS = CENTER + 0.5*SIDES.*[-1 -1 ; 1 -1 ; 1 1 ; -1 1] ;
    d = pkg.geometry.distance.toPolygon(P,PTS) ;

end

