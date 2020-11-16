classdef Rectangle < pkg.geometry.levelset.LevelSet
% Level set representing a Rectangle

methods
    function this = Rectangle(box,sz)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        if nargin<1 ; return ; end
        center = mean(box,1) ;
        if nargin<2 ; sz = range(box,1) ; end
        this.Function = @(p)pkg.geometry.distance.toRectangle(p,center,sz) ;
        this.Kinks = center + sz/2.*[-1 -1 ; 1 -1 ; 1 1 ; -1 1] ;
        this.BoundingBox = this.Kinks([1,3],:) ;
        this.EdgeFcns = pkg.geometry.levelset.polylineEdgeFunctions(this.Kinks([1:end,1],:)) ;
    end
end
    
end

