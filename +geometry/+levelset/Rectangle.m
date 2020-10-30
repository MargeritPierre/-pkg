classdef Rectangle < pkg.geometry.levelset.LevelSet
% Level set representing a Rectangle

methods
    function this = Rectangle(center,sides)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = @(p)pkg.geometry.distance.toRectangle(p,center,sides) ;
        this.Kinks = center + sides/2.*[-1 -1 ; 1 -1 ; 1 1 ; -1 1] ;
        this.BoundingBox = this.Kinks([1,3],:) ;
        this.EdgeFcns = pkg.geometry.levelset.polylineEdgeFunctions(this.Kinks([1:end,1],:)) ;
    end
end
    
end

