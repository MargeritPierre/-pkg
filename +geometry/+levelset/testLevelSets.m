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
    lvlst = ls1 - ls2 - ls3 + ls4 - ls5 + ls6 ;
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
    dx = norm(range(lvlst.BoundingBox,1))/50 ;
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
        






    



