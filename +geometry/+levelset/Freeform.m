classdef Freeform < pkg.geometry.levelset.LevelSet
% Level set representing a Freeform surface/volume

methods
    function this = Freeform(fcn,bbox)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = fcn ;
        this.BoundingBox = bbox ;
        cutbox = pkg.geometry.levelset.Box(bbox) ;
        %[this.EdgeFcns,this.Kinks] = this.intersectEdges(cutbox) ;
    % Cut with the boundingbox and copy the content
        cutlvlst = this & cutbox ;
        this.Function = cutlvlst.Function ;
        this.EdgeFcns = cutlvlst.EdgeFcns ;
        this.Kinks = cutlvlst.Kinks ;
%         this.cleanContour() ;
    end
end
    
end

