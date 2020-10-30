classdef Polyline < pkg.geometry.levelset.LevelSet
% Level set representing a Polyline

methods
    function this = Polyline(points)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = @(p)pkg.geometry.distance.toPolyline(p,points) ;
        this.BoundingBox = [min(points,[],1) ; max(points,[],1)] ;
        this.Kinks = points ;
        this.EdgeFcns = pkg.geometry.levelset.polylineEdgeFunctions(points) ;
    end
end
    
end

