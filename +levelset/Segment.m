classdef Segment < pkg.levelset.LevelSet
% Level set representing a Segment


    
    methods
        function this = Segment(p1,p2)
        % Class construtor
            this = this@pkg.levelset.LevelSet() ;
            this.Function = @(p)pkg.levelset.LevelSet.dsegment(p,p1,p2) ;
            this.BoundingBox = [min(p1,p2) ; max(p1,p2)] ;
            this.Kinks = [p1 ; p2] ;
            this.EdgeFcns{1} = @(t)p1.*(1-t(:)) + p2.*t(:) ;
        end
    end
end

