function d = toPolygon(P,PTS)
%TOPOLYGON SIGNED Distance from a list of points P to a CLOSED polygon PG = [PTS(1)->PTS(2)->...]
% inputs:
%   - P [nP nCoord]
%   - PTS [nPinPolygon nCoord 1]
% output: 
%   - d [nP 1 1] signed distance (d<0 inside, d>0 outside)

    if isempty(P) ; d = [] ; return ; end

    % Distance to a polyline
    d = pkg.geometry.distance.point.toPolyline(P,PTS([1:end,1],:)) ;
    % Signed distance
    d = d.*(0.5-inpolygon(P(:,1),P(:,2),PTS(:,1),PTS(:,2)))*2 ;

end

