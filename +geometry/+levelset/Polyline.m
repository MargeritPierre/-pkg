classdef Polyline < pkg.geometry.levelset.LevelSet
% Level set representing a Polyline

methods
    function this = Polyline(points,edgFcn,t)
    % Class construtor
        if nargin<2 || isempty(edgFcn) ; edgFcn = 'segments' ; end
        if nargin<3 || isempty(t) ; t = 0 ; end
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = @(p)pkg.geometry.distance.point.toPolyline(p,points)-t/2 ;
        this.BoundingBox = [min(points,[],1) ; max(points,[],1)] ;
        this.Kinks = points ;
        this.EdgeFcns = pkg.geometry.levelset.polylineEdgeFunctions(points,edgFcn) ;
    end
end
    
end

