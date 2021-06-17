classdef Cylinder < pkg.geometry.levelset.LevelSet
% Level set representing a Cylinder

methods
    function this = Cylinder(radius,pts,isfinite)
    % Class construtor
        if nargin<2 ; pts = [0 0 0;0 0 1] ; end
        if nargin<3 ; isfinite = false ; end
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = @(p)pkg.geometry.distance.toCylinder(p,pts,isfinite) ;
        this.Kinks = [] ;
        if isfinite
            this.BoundingBox = sort(pts,1) + radius*[-1 -1 -1 ; 1 1 1] ; % approximated!
            this.EdgeFcns{1} = @(t)pts(1,:) + radius*[cos(2*pi*t(:)) sin(2*pi*t(:)) 0*t(:)] ;
            this.EdgeFcns{2} = @(t)pts(2,:) + radius*[cos(2*pi*t(:)) sin(2*pi*t(:)) 0*t(:)] ;
        end
    end
end
    
end

