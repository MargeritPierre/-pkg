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
    % Build the boundaries
        this = this.buildEdges ;
    % Set the distance function
        this.Function = @this.distanceFun ;
    end
    
    function d = distanceFun(this,p)
    % Signed distance function
    % Comply with a custom bounding box 
    % !!!! distance in PIXELS, the gradient is NOT RIGHT !!!
        %p = p-this.BoundingBox(1,:) ;
        %p = p./range(this.BoundingBox,1) ;
    % Image indices (i,j) = (y,x)
        IDX = flip(p,2) ; %.*size(this.LevelSetImage) ;
    % Interpolation
        ORDER = 1 ; % linear interpolation
        EXTRAP = true ; % extrapolate ouside the image
        d = pkg.data.interp(this.LevelSetImage,IDX,ORDER,EXTRAP) ;
    end
    
    function this = buildEdges(this)
    % Build the levelset edges using image contour
    % (like a marching square algo)
    % First find where the levelset crosses 0
        img = this.LevelSetImage ;
        si = sign(img) ;
        si = cat(3 , si(1:end-1,1:end-1) ...
                    ,si(2:end,1:end-1) ...
                    ,si(1:end-1,2:end) ...
                    ,si(2:end,2:end)) ;
        cross = any(diff(si,1,3),3) ;
        nElems = sum(cross(:)) ;
    % Create a quad mesh from the crossing locations
        % crossing locations
        [ii,jj] = find(cross) ;
        % build the quad corners
        p1 = [jj(:) ii(:)] ;
        p2 = p1 + [1 0] ; p3 = p2 + [0 1] ; p4 = p1 + [0 1] ;
        X = [p1 ; p2 ; p3 ; p4] ;
        % unique nodes
        [X,~,ic] = unique(X,'rows') ;
        % add the levelset in the third coordinate
        X = [X img(sub2ind(size(img),X(:,2),X(:,1)))] ;
        % quad elems
        idx = (1:nElems)' + (0:3)*nElems ;
        idx = reshape(ic(idx),size(idx)) ;
        % mesh
        elems = pkg.geometry.mesh.elements.ElementTable(...
                    'Types',pkg.geometry.mesh.elements.base.Quadrangle ...
                    ,'Indices',padarray(idx,[0 1],1,'pre')) ;
        mesh = pkg.geometry.mesh.Mesh('Nodes',X,'Elems',elems) ;
    % Slice the mesh at x(:,3)==0
        isoMesh = mesh.cut(@(x)x(:,3)).ON ;
        isoMesh.nCoord = 2 ; % constrain in the plane
    % Separate the mesh curves
        [~,crvs] = isoMesh.boundaryCurves ;
    % Convert to edge functions
        this.EdgeFcns = {} ;
        for ee = 1:numel(crvs)
            points = isoMesh.Nodes(crvs{ee},:) ;
            this.EdgeFcns{end+1} = pkg.geometry.levelset.polylineEdgeFunctions(points,'full') ;
        end
    end
end
    
end

