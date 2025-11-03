classdef Sphere < pkg.geometry.levelset.LevelSet
% Level set representing a Sphere

methods
    function this = Sphere(center,radius)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        center(end+1:3) = 0 ; % force 3D coordinates
        this.Function = @(p)pkg.geometry.distance.point.toSphere(p,center,radius) ;
        this.BoundingBox = center + radius.*[-1 ; 1] ;
        this.Kinks = [] ;
        this.EdgeFcns = {} ;
        this.SurfMeshFcn = @(dx)pkg.geometry.mesh.shapes.sphere(dx/radius).scale(radius).move(center) ;
    end
end
    
end

