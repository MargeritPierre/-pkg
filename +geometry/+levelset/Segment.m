classdef Segment < pkg.geometry.levelset.LevelSet
% Level set representing a Segment

methods
    function this = Segment(p1,p2,t)
    % Class construtor
        if nargin<3 || isempty(t) ; t = 0 ; end
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = @(p)pkg.geometry.distance.point.toSegment(p,p1,p2)-t/2 ;
        this.BoundingBox = [min(p1,p2) ; max(p1,p2)] ;
        this.Kinks = [p1 ; p2] ;
        this.EdgeFcns{1} = @(t)p1.*(1-t(:)) + p2.*t(:) ;
    end
end
    
end

