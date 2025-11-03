classdef Plane < pkg.geometry.levelset.LevelSet
% Level set representing an infinite plane

methods
    function this = Plane(origin,normal)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        origin(end+1:3) = 0 ; % force 3D coordinates
        normal(end+1:3) = 0 ; % force 3D coordinates
        this.Function = @(p)pkg.geometry.distance.point.toPlane(p,origin,normal) ;
        this.BoundingBox = origin + inf.*[-1 ; 1] ;
        this.Kinks = [] ;
        this.EdgeFcns = {} ;
    end
end
    
end

