classdef Ellipse < pkg.geometry.levelset.LevelSet
% Level set representing an Ellipse

methods
    function this = Ellipse(center,semiaxes,rotation)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        this.Function = @(p)pkg.geometry.distance.point.toEllipse(p,center,semiaxes,rotation) ;
        this.Kinks = [] ;
        this.EdgeFcns{1} = @(t)center + [semiaxes(1)*cos(2*pi*t(:)) semiaxes(2)*sin(2*pi*t(:))]*[cos(rotation) sin(rotation) ; -sin(rotation) cos(rotation)] ;
        t = [-atan(semiaxes(2)/semiaxes(1)*tan(rotation)) ; atan(semiaxes(2)/semiaxes(1)/tan(rotation))] ;
        t = [t ; t+pi]/2/pi ;
        extrPts = this.EdgeFcns{1}(t) ;
        this.BoundingBox = [min(extrPts,[],1) ; max(extrPts,[],1)] ;
    end
end
    
end

