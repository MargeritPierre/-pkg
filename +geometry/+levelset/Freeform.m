classdef Freeform < pkg.geometry.levelset.LevelSet
% Level set representing a Freeform surface/volume

methods
    function this = Freeform(fcn,bbox,surface)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = fcn ;
        this.BoundingBox = bbox ;
    % Cut with the boundingbox
        cutbox = pkg.geometry.levelset.Box(bbox) ;
        %[this.EdgeFcns,this.Kinks] = this.intersectEdges(cutbox) ;
        cutlvlst = this & cutbox ;
        this.Kinks = cutlvlst.Kinks ;
        this.EdgeFcns = cutlvlst.EdgeFcns ;
        if nargin>=3 && surface
            this.Function = @(p)abs(fcn(p)) ;
            this = this.cleanContour() ;
            this.Function = fcn ;
        else
            this.Function = cutlvlst.Function ;
        end
    end
end
    
end

