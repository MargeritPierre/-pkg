classdef Circle < pkg.geometry.levelset.LevelSet
% Level set representing a Circle

methods
    function this = Circle(center,radius)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = @(p)pkg.geometry.distance.toCircle(p,center,radius) ;
        this.BoundingBox = center + radius*[-1 -1 ; 1 1] ;
        this.Kinks = [] ;
        this.EdgeFcns{1} = @(t)center + radius*[cos(2*pi*t(:)) sin(2*pi*t(:))] ;
    end
end
    
end

