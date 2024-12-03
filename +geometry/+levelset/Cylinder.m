classdef Cylinder < pkg.geometry.levelset.LevelSet
% Level set representing a Cylinder

methods
    function this = Cylinder(radius,pts,isfinite)
    % Class construtor
        if nargin<2 ; pts = [0 0 0;0 0 1] ; end
        if nargin<3 ; isfinite = true ; end
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = @(p)pkg.geometry.distance.point.toCylinder(p,radius,pts,~isfinite) ;
        this.Kinks = [] ;
        if isfinite
            this.BoundingBox = sort(pts,1) + radius*[-1 -1 -1 ; 1 1 1] ; % approximated!
            axvec = diff(pts,1,1) ; L = norm(axvec) ;
            uaxvec = axvec./L ; % unit cylinder axis vector
            rotax = pkg.math.vectprod([0 0 1],uaxvec,2) ;
            ROT = pkg.math.rotmat(-rotax) ;
            this.EdgeFcns{1} = @(t)pts(1,:) + radius*[cos(2*pi*t(:)) sin(2*pi*t(:)) 0*t(:)]*ROT ;
            this.EdgeFcns{2} = @(t)pts(2,:) + radius*[cos(2*pi*t(:)) sin(2*pi*t(:)) 0*t(:)]*ROT ;
            this.SurfMeshFcn = @(dx)pkg.geometry.mesh.shapes.cylinder(L/radius,dx).scale(radius).applyTransform(ROT).move(pts(1,:)) ;
        else 
            this.BoundingBox = Inf(1,3).*[-1;1] ;
        end
    end
end
    
end

