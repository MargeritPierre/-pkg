%% UNIT TESTS FOR THE LEVELSET CLASS

%% BASE SHAPES

% Different base shapes
    lvlst = pkg.geometry.levelset.Point([0 1]) ; % (position)
    lvlst = pkg.geometry.levelset.Line([-1 1],[1 -1]) ; % (p1,p2)
    lvlst = pkg.geometry.levelset.Segment([0 0],[1 1]) ; % (p1,p2)
    lvlst = pkg.geometry.levelset.Polyline([-1 1 ; 0 -1 ; 1 1 ; 0 0]) ; % (points)
    lvlst = pkg.geometry.levelset.Polygon([-1 1 ; 0 -1 ; 1 1 ; 0 0]) ; % (points)
    lvlst = pkg.geometry.levelset.Rectangle([1 1],[2 3]) ; % (center,sides)
    lvlst = pkg.geometry.levelset.Circle([0,1],1) ; % (center,radius)
    lvlst = pkg.geometry.levelset.Ellipse([0,1],[1 3],45/180*pi) ; % (center,semiaxes,rotation)

% Enhance BBox
    if ~isempty(lvlst.BoundingBox) ; lvlst.BoundingBox = lvlst.BoundingBox + [-1 -1 ; 1 1]*1 ; end

% Plot the level set
    clf
    axis equal
    axis tight
    tic ; h = plot(lvlst) ; toc
    colorbar
    
    
%% IMAGE MASK

% Round mask
    R = 500 ; r = R*3/4 ; MASK = (-R:R)'.^2 + (-R:R).^2 <= r.^2 ; 
    tic ; lvlst = pkg.geometry.levelset.Image(MASK) ; toc
    
% Enhance BBox
    if ~isempty(lvlst.BoundingBox) ; lvlst.BoundingBox = lvlst.BoundingBox + [-1 -1 ; 1 1]*0 ; end
    
% Plot the level set
    clf
    axis equal
    axis tight
    tic ; h = plot(lvlst) ; toc
    colorbar
    
    
%% BOOLEAN OPERATIONS

% Create base shapes
    ls1 = pkg.geometry.levelset.Rectangle([0 0],[4 4]) ;
    ls2 = pkg.geometry.levelset.Ellipse([2 0],[1 2],40/180*pi) ;
    ls3 = pkg.geometry.levelset.Ellipse([-2 -1],[1 2],-40/180*pi) ;
    ls4 = pkg.geometry.levelset.Circle([0 0],1.5) ;
    ls5 = pkg.geometry.levelset.Rectangle([0 0],[1 1]) ;
    ls6 = pkg.geometry.levelset.Polygon([-1 1 ; 0 -1 ; 1 1 ; 0 0]) ;
    
% Assemble them
    profile on
    tic
    lvlst = ls1 - ls2 - ls3 + ls4 - ls5 ; + ls6 ;
    toc
    profile off
    
% Enhance BBox
    if ~isempty(lvlst.BoundingBox) ; lvlst.BoundingBox = lvlst.BoundingBox + [-1 -1 ; 1 1]*0.1 ; end

% Plot the level set
    clf
    axis equal
    axis tight
    tic ; h = plot(lvlst) ; toc
    %caxis([-1 1]*0.01)
    
%% MIX SHAPES & IMAGE

% Round mask
    R = 1000 ; r = R*3/4 ; MASK = (-R:R)'.^2 + (-R:R).^2 <= r.^2 ; 
    lvlst = pkg.geometry.levelset.Image(MASK) ;
    
% Shapes
    lvlst = lvlst - pkg.geometry.levelset.Rectangle(R*[1.5 1.5],R*[1 1]) ;
    lvlst = lvlst - pkg.geometry.levelset.Line(R*[.5 0],R*[.5 2]) ;
    
% Enhance BBox
    if ~isempty(lvlst.BoundingBox) ; lvlst.BoundingBox = lvlst.BoundingBox + [-1 -1 ; 1 1]*0.1 ; end

% Plot the level set
    clf
    axis equal
    axis tight
    tic ; h = plot(lvlst) ; toc
    %caxis([-1 1]*0.01)
    
    
%% LEVEL SET CONTOUR DISCRETIZATION
    dx = norm(range(lvlst.BoundingBox,1))/50 ;
    dx = @(p)dx-1.5*(p(:,2)/norm(range(lvlst.BoundingBox,1))).^2 ;
    P = lvlst.discretizeContour(dx) ;
    tag = 'contour' ;
    delete(findobj(gca,'tag',tag)) ;
    c = plot(P(:,1),P(:,2),'.k','markersize',15,'tag',tag) ;
    
%% LEVELSET SKELETON
    dx = norm(range(lvlst.BoundingBox,1))/500 ;
    skel = lvlst.skeleton(dx) ;
    ps = plot(skel) ;

%% POPULATE THE LEVELSET
    dx = 0.01 ;
    tag = 'population' ;
    delete(findobj(gca,'tag',tag)) ;
    P = lvlst.populate(dx,'iso',@(P)dx*(1 + abs(P(:,2))*2)) ;
    %P = [P ; lvlst.discretizeContour(dx)] ;
    c = plot(P(:,1),P(:,2),'.k','markersize',8,'tag',tag) ;
    
    
%% BUILD A MESH
    delete(findobj(gca,'type','hggroup'))
    
    profile on
    mesh = pkg.geometry.mesh.distMesh(lvlst ...
                                        ... ,'h0',3 ...
                                        ... ,'debug',false ...
                                        ,'showMesh',true ...
                                        ... ,'tooLongThrs',1.4 ...
                                        ... ,'tooShortThrs',0.6 ...
                                        ... ,'maxCount',20 ...
                                        ... ,'deltat',1 ...
                                        ... ,'qmin',0.3 ...
                                        ... ,'Fscale',1.2 ...
                                        ... ,'bndCons',false ...
                                        ... ,'t_dmax',0 ...
                                        ) ;
    profile off
    
    drawnow ;
    delete(findobj(gca,'type','hggroup'))
    pl = plot(mesh) ;
    

%% MESH A POLYGON with variable segment length

    N = 50 ; t = linspace(0,1-1/(2*N),N).^1.4' ; t = [-flip(t(:)) ; t(:)] ;
    P = [1.0*cos(pi*t(:)) 1.0*sin(pi*t(:))] ;
    P = unique(P,'rows','stable') ;
    lvlst = pkg.geometry.levelset.Polygon(P) ;
    dens = pkg.geometry.density.Polygon(P) ;
    h0 = min(dens.Le) ;
    fh = @dens.evalAt ;
    
    clf
    axis equal
    axis tight
    h = plot(lvlst) ;
    
    mesh = pkg.geometry.mesh.distMesh(lvlst ...
                        ,'h0',h0 ...
                        ,'fh',fh ...
                        ,'pfix',P ...
                        ,'showMesh',true) ;
    
    drawnow ;
        delete(findobj(gca,'type','hggroup'))
        delete(h)
    pl = plot(mesh) ;
        pl.FaceColor = 'flat' ;
        pl.Faces.FaceVertexCData = mesh.elemSize ;
        

%% ADVANCING FRONT MESH
    tag = 'advancingfront' ;
    delete(findobj(gca,'tag',tag)) ;
    mesh = pkg.geometry.mesh.Mesh() ;
    pl = plot(mesh,'Tag',tag) ;
    pl.VisibleNodes = 'all' ;
    pl.Nodes.MarkerSize = 10 ;
    
    fd = @lvlst.Function ;
    fgrad = @lvlst.gradient ;
    
    h0 = norm(range(lvlst.BoundingBox,1))/50 ;
    fh = @(p)repmat(h0,[size(p,1) 1]) ;
    
    % Start from the discretized contour
        mesh.Nodes = lvlst.discretizeContour(fh);
        lastNodes = 1:mesh.nNodes ;
    
    % Edges on the front: initialize
        mesh.Elems.Types = pkg.geometry.mesh.elements.base.Triangle ;
%         mesh.Elems.NodeIdx = delaunay(mesh.Nodes) ;
%         mesh.Elems = mesh.Elems.subpart(lvlst.Function(mesh.centroid)<0) ;
    
    maxIt = 2*norm(range(lvlst.BoundingBox,1))/h0 ;
    deps = 1e-3 ;
    it = 0 ;
    while ~isempty(lastNodes) && it<maxIt
        it = it+1 ;
        p = mesh.Nodes(lastNodes,:) ;
    % Keep a copy of the current points
        p0 = p ;
        d0 = fd(p) ;
    % New layer of points
        for iii = 1:10
        % Direction of the gradient
            grad = fgrad(p) ;
        % Target edge length
            ht = fh((p0+p)/2) ;
        % Current edge length
            h = d0-fd(p) ; %sqrt(sum((p0-p).^2,2)) ;
        % Move accordingly
            p = p-grad.*(ht-h) ;
        end
    % Cull points for which the signed distance increased
        d = fd(p) ;
        p(d>d0,:) = [] ;
        d0 = d(d>d0) ;
    % Add Nodes to the mesh
        lastNodes = mesh.nNodes + (1:size(p,1)) ;
        mesh.Nodes = [mesh.Nodes ; p] ;
    % Display
        if 1 ; pl.update ; drawnow ; end
    end
    
    pl.update
    
    





    



