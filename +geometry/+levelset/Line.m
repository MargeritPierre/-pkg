classdef Line < pkg.geometry.levelset.LevelSet
% Level set representing a Line
    
methods
    function this = Line(p1,p2,t)
    % Class construtor
        if nargin<3 || isempty(t) ; t = 0 ; end
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = @(p)pkg.geometry.distance.point.toLine(p,p1,p2)-t/2 ;
        this.BoundingBox = [] ;
        this.Kinks = [] ;
        this.EdgeFcns{1} = @(t)p1.*(1-t(:)) + p2.*t(:) ;
    end
end

end

