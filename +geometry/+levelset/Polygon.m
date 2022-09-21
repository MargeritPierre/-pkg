classdef Polygon < pkg.geometry.levelset.LevelSet
% Level set representing a Polygon

methods
    function this = Polygon(points,edgFcn)
    % Class construtor
        if nargin<2 ; edgFcn = 'segments' ; end
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = @(p)pkg.geometry.distance.point.toPolygon(p,points) ;
        this.BoundingBox = [min(points,[],1) ; max(points,[],1)] ;
        this.Kinks = points ;
        this.EdgeFcns = pkg.geometry.levelset.polylineEdgeFunctions(points([1:end,1],:),edgFcn) ;
    end
end
    
end

