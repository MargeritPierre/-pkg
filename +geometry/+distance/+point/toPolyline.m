function d = toPolyline(P,PTS)
%TOPOLYLINE Distance from a list of points P to a polyline PL = [PTS(1)->PTS(2)->...]
% The polyline MAY be closed (if PTS(1)==PTS(end))
% inputs:
%   - P [nP nCoord]
%   - PTS [nPinPolyline nCoord 1]
% output: 
%   - distance d [nP 1 1] 

    if isempty(P) ; d = [] ; return ; end

    PTS = permute(PTS,[3 2 1]) ;
    dsegments = pkg.geometry.distance.point.toSegment(P,PTS(:,:,1:end-1),PTS(:,:,2:end)) ;
    d = min(dsegments(:,:),[],2) ;

end

