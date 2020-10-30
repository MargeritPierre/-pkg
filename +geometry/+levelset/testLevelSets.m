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
    lvlst = pkg.geometry.levelset.Image(MASK) ;
    
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
    %P = discretizeContour(lvlst) ;
    P = lvlst.discretizeContour(@(p)dx*(1-0.2*p(:,1))) ;
    c = plot(P(:,1),P(:,2),'.k','markersize',15) ;
    
%% LEVELSET SKELETON
    skel = lvlst.skeleton ;
    ps = plot(skel) ;

%% POPULATE THE LEVELSET
    dx = 0.011 ;
    tag = 'population' ;
    delete(findobj(gca,'tag',tag)) ;
    P = lvlst.populate(dx,'iso',@(P)dx*(1 + abs(P(:,2))*2)) ;
    %P = [P ; lvlst.discretizeContour(dx)] ;
    c = plot(P(:,1),P(:,2),'.k','markersize',8,'tag',tag) ;
    
    
%% BUILD A MESH
    delete(findobj(gca,'type','hggroup'))
    
    %profile on
    mesh = pkg.geometry.mesh.distMesh(lvlst) ;
    %profile off
    
    drawnow ;
    delete(findobj(gca,'type','hggroup'))
    pl = plot(mesh) ;
    



