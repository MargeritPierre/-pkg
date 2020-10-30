classdef Line < pkg.geometry.levelset.LevelSet
% Level set representing a Line
    
methods
    function this = Line(p1,p2)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = @(p)pkg.geometry.distance.toLine(p,p1,p2) ;
        this.BoundingBox = [] ;
        this.Kinks = [] ;
        this.EdgeFcns{1} = @(t)p1.*(1-t(:)) + p2.*t(:) ;
    end
end

end

