classdef Point < pkg.geometry.levelset.LevelSet
% Level set representing a point
    
methods
    function this = Point(pos)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = @(p)this.dpoint(p,pos) ;
        this.BoundingBox = [] ;
        this.Kinks = pos ;
        this.EdgeFcns{1} = @(t)repmat(pos,size(t(:))) ;
    end
end
    
end

