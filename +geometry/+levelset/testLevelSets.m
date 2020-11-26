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
                                        ... ,'h0',10 ...
                                        ... ,'pfix',lvlst.discretizeContour ...
                                        ... ,'p0',mesh.Nodes ...
                                        ... ,'debug',true ...
                                         ,'showMesh',true ...
                                        ... ,'tooLongThrs',1.4 ...
                                        ... ,'tooShortThrs',0.6 ...
                                        ... ,'maxCount',20 ...
                                        ... ,'deltat',1 ...
                                        ... ,'qmin',0.1 ...
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
        
        
%% UPDATE A MESH WITH A CHANGING DENSITY  

% Parameters
    hmax = 0.5 ;
    hmin = 0.01 ;
    l = 1.5 ;
    pts = [0 1 0 ; 1 -3 6] ; % [time x y]
    t = linspace(0,1,100)*range(pts(:,1)) + min(pts(:,1)) ;
    constArgs = { ...
             'h0',hmin ...
             ,'qmin',0*0.1 ...
             ,'bndCons',true ...
             ...,'tooShortThrs',0 ...
             ...,'tooLongThrs',Inf ...
             ,'showMesh',false ...
            } ;
% Init display
    delete(findobj(gca,'type','hggroup'))
    pl = pkg.geometry.mesh.MeshPlot('Parent',gca) ;
% Loop
profile on
    for tt = 1:numel(t)
    % Interpolate the point
        pt = interp1(pts(:,1),pts(:,2:3),t(tt)) ;
    % Create the density function
        dens = pkg.geometry.density.Point(pt,[0 hmin ; l hmax]) ;
        fh = @dens.evalAt ;
    % Update the mesh points
        if tt==1 % Init the mesh if needed
            p0 = lvlst.populate(hmin,'iso',fh) ;
            mesh = pkg.geometry.mesh.Mesh('Nodes',p0) ;
        else
            pfix = lvlst.discretizeContour(fh) ;
            p0 = mesh.Nodes(~mesh.boundaryNodes,:) ; %lvlst.populate(hmin,'iso',fh) ;
        end
    % Relaxation
        tic
        mesh = pkg.geometry.mesh.distMesh(lvlst ...
                                            ,'fh',fh ...
                                            ,'p0',p0 ...
                                            ,'pfix',pfix ...
                                            ,constArgs{:} ...
                                            ) ;
        toc
    % Display
        pl.Mesh = mesh ;
        drawnow ;
    end
profile off


        






    



