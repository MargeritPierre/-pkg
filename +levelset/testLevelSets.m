%% UNIT TESTS FOR THE LEVELSET CLASS

%% BASE SHAPES

% Different base shapes
    lvlst = pkg.levelset.Point([0 1]) ; % (position)
    lvlst = pkg.levelset.Line([-1 1],[1 -1]) ; % (p1,p2)
    lvlst = pkg.levelset.Segment([0 0],[1 1]) ; % (p1,p2)
    %lvlst = pkg.levelset.LevelSet.Polyline([-1 1 ; 0 -1 ; 1 1 ; 0 0]) ; % (points)
    %lvlst = pkg.levelset.LevelSet.Polygon([-1 1 ; 0 -1 ; 1 1 ; 0 0]) ; % (points)
    %lvlst = pkg.levelset.LevelSet.Rectangle([1 1],[2 3]) ; % (center,sides)
    %lvlst = pkg.levelset.LevelSet.Circle([0,1],1) ; % (center,radius)
    %lvlst = pkg.levelset.LevelSet.Ellipse([0,1],[1 3],45/180*pi) ; % (center,semiaxes,rotation)

% Enhance BBox
    if ~isempty(lvlst.BoundingBox) ; lvlst.BoundingBox = lvlst.BoundingBox + [-1 -1 ; 1 1]*1 ; end

% Plot the level set
    clf
    axis equal
    axis tight
    tic ; h = plot(lvlst) ; toc
    
    
%% BOOLEAN OPERATIONS

% Create base shapes
    ls1 = pkg.levelset.LevelSet.Rectangle([0 0],[4 4]) ;
    ls2 = pkg.levelset.LevelSet.Ellipse([2 0],[1 2],40/180*pi) ;
    ls3 = pkg.levelset.LevelSet.Ellipse([-2 -1],[1 2],-40/180*pi) ;
    ls4 = pkg.levelset.LevelSet.Circle([0 0],1.5) ;
    ls5 = pkg.levelset.LevelSet.Rectangle([0 0],[1 1]) ;
    
% Assemble them
    %profile on
    tic
    lvlst = ls1 - ls2 - ls3 + ls4 - ls5 ;
    toc
    %profile viewer
    
% Enhance BBox
    if ~isempty(lvlst.BoundingBox) ; lvlst.BoundingBox = lvlst.BoundingBox + [-1 -1 ; 1 1]*1 ; end

% Plot the level set
    clf
    axis equal
    axis tight
    tic ; h = plot(lvlst,500) ; toc
    %caxis([-1 1]*0.01)
    
% Discretize contour
    P = discretizeContour(lvlst) ;
    c = plot(P(:,1),P(:,2),'.k','markersize',15) ;
    
%% LEVEL SET CONTOUR

    [C,ch] = contour(h(1).XData,h(1).YData,h(1).ZData,[1 1]*0,'k') ;




