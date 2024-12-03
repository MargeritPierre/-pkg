classdef Box < pkg.geometry.levelset.LevelSet
% Level set representing a Box

methods
    function this = Box(box,sz)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        if nargin<1 ; return ; end
        % Force 3D coordinates
        box(:,end+1:3) = 0 ; 
        % Process input
        center = mean(box,1) ;
        if nargin<2 ; sz = range(box,1) ; end
        sz(end+1:3) = 0 ; 
        % Build
        this.Function = @(p)pkg.geometry.distance.point.toBox(p,center,sz) ;
        this.Kinks = center + sz/2.*[-1 -1 -1 ; 1 -1 -1 ; 1 1 -1 ; -1 1 -1 ; -1 -1 1 ; 1 -1 1 ; 1 1 1 ; -1 1 1] ;
        this.BoundingBox = center + sz.*[-1;1]/2 ;
        edg = [1 2 ; 2 3 ; 3 4 ; 4 1 ; 5 6 ; 6 7 ; 7 8 ; 8 5 ; 1 5 ; 2 6 ; 3 7 ; 4 8] ;
        this.EdgeFcns = cellfun(@(ee)pkg.geometry.levelset.polylineEdgeFunctions(this.Kinks(ee,:)),num2cell(edg,2)) ;
        this.SurfMeshFcn = @(dx)pkg.geometry.mesh.shapes.box(sz,dx).move(center-sz/2) ;
    end
end
    
end

