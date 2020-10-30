classdef Image < pkg.geometry.levelset.LevelSet
% Level set representing a BINARY Image
%   where white pixels (1's) sample the geometry

properties (SetAccess = immutable)
    LevelSetImage
end

methods
    function this = Image(IMG,bbox)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
    % Set common levelset properties
        if nargin<2 ; bbox = [1 1 ; size(IMG)] ; end
        this.BoundingBox = bbox ;
        this.Kinks = [] ;
        this.EdgeFcns = {} ;
    % Build the levelset image
        ROI = logical(IMG) ;
        DIST = bwdist(ROI)-bwdist(~ROI) ;
        this.LevelSetImage = double(DIST + ROI.*0.5 - ~ROI*0.5) ;
    % Set the distance function
        this.Function = @this.distanceFun ;
    end
    
    function d = distanceFun(this,p)
    % Signed distance function
    % Comply with a custom bounding box 
    % !!!! distance in PIXELS, the gradient is NOT RIGHT !!!
        p = p-this.BoundingBox(1,:) ;
        p = p./range(this.BoundingBox,1) ;
    % Image indices (i,j) = (y,x)
        IDX = flip(p,2).*size(this.LevelSetImage) ;
    % Interpolation
        ORDER = 1 ; % linear interpolation
        EXTRAP = true ; % extrapolate ouside the image
        d = pkg.data.interp(this.LevelSetImage,IDX,ORDER,EXTRAP) ;
    end
end
    
end

