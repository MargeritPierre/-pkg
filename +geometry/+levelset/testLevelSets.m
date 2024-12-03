%% UNIT TESTS FOR THE LEVELSET CLASS

%% BASE SHAPES

% 2D base shapes
    lvlst = pkg.geometry.levelset.Point([0 1]) ; % (position)
    lvlst = pkg.geometry.levelset.Line([-1 1],[1 -1]) ; % (p1,p2)
    lvlst = pkg.geometry.levelset.Segment([0 0],[1 1]) ; % (p1,p2)
    lvlst = pkg.geometry.levelset.Polyline([-1 1 ; 0 -1 ; 1 1 ; 0 0]) ; % (points)
    lvlst = pkg.geometry.levelset.Polygon([-1 1 ; 0 -1 ; 1 1 ; 0 0]) ; % (points)
    lvlst = pkg.geometry.levelset.Rectangle([1 1],[2 3]) ; % (center,sides)
    lvlst = pkg.geometry.levelset.Circle([0,1],1) ; % (center,radius)
    lvlst = pkg.geometry.levelset.Ellipse([0,1],[1 3],45/180*pi) ; % (center,semiaxes,rotation)
 
% 3D base shapes
    lvlst = pkg.geometry.levelset.Sphere([0,1,0],1) ; % (center,radius)
    lvlst = pkg.geometry.levelset.Box([0,1,1.5],[1 2 1.5]) ; % (center,sides)
    %lvlst = pkg.geometry.levelset.Cylinder(.74,[1.3 2.2 0 ; 0 0 3.2],true) ; % (radius,axis end pts,is finite)
    
% Enhance BBox
    if ~isempty(lvlst.BoundingBox) ; lvlst.BoundingBox = lvlst.BoundingBox + [-1 ; 1]*1 ; end

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
    
%% FREEFORM

    fcn = @(x)sin(x(:,1)).*cos(x(:,2)) + sin(x(:,2)).*cos(x(:,3)) + sin(x(:,3)).*cos(x(:,1)) ;
    bbox = pi.*[-1 ; 1].*[1 1 1] ;
    lvlst = pkg.geometry.levelset.Freeform(fcn,bbox)
    t = .25 ; lvlst = pkg.geometry.levelset.Freeform(@(p)abs(fcn(p))-t,bbox) ;
    
% Enhance BBox
    %if ~isempty(lvlst.BoundingBox) ; lvlst.BoundingBox = lvlst.BoundingBox + [-1 ; 1]*1 ; end

% Plot the level set
    clf
    axis equal
    axis tight
    tic ; h = plot(lvlst) ; toc
    colorbar
    
%% 2D BOOLEAN OPERATIONS

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
    
    
%% 3D BOOLEAN OPERATIONS

% Create base shapes
    ls1 = pkg.geometry.levelset.Sphere([0,1,0],1) ; % (center,radius)
    ls2 = pkg.geometry.levelset.Sphere([0,0,0],1) ; % (center,radius)
    ls3 = pkg.geometry.levelset.Cylinder(.4,[0 2 0 ; 0 -1 0],true) ; % (radius,axis end pts,is finite)
    ls4 = pkg.geometry.levelset.Box([0,.5,0],[3 1 1]) ; % (center,sides)
    
% Assemble them
    profile on
    tic
    lvlst = (ls1 + ls2) - ls4 ;
    toc
    profile off
    
% Enhance BBox
    if ~isempty(lvlst.BoundingBox) ; lvlst.BoundingBox = lvlst.BoundingBox + [-1 ; 1]*0.1 ; end

% Plot the level set
    clf
    axis equal
    axis tight
    tic ; h = plot(lvlst) ; toc
    %caxis([-1 1]*0.01)
    
%% SMOOTH INTERSECTIONS
% see https://iquilezles.org/articles/smin/

% 2D example
ls1 = pkg.geometry.levelset.Rectangle([0 0],[1,1]) ; % (center,sides)
ls2 = pkg.geometry.levelset.Circle(7/10*[1 1],.5) ; % (center,radius)

% 3D example
ls1 = pkg.geometry.levelset.Cylinder(.5,[-2 0 0 ; 2 0 0]) ; % (radius,endpoints)
ls2 = pkg.geometry.levelset.Cylinder(.5,[0 -2 .5 ; 0 2 .5]) ; % (center,radius)

lvlst = merge(ls1,ls2,1/10) ;
    
% Enhance BBox
    if ~isempty(lvlst.BoundingBox) ; lvlst.BoundingBox = lvlst.BoundingBox + [-1 ; 1]*0.1 ; end

% Plot the level set
    clf
    axis equal
    axis tight
    tic ; h = plot(lvlst) ; toc
    caxis([-1 1]*0.01)


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
    %dx = @(p)dx-1.5*(p(:,2)/norm(range(lvlst.BoundingBox,1))).^2 ;
    P = lvlst.discretizeContour(dx) ;
    tag = 'contour' ;
    delete(findobj(gca,'tag',tag)) ;
    P(:,end+1:3) = 0 ; % 3D coordinates
    c = plot3(P(:,1),P(:,2),P(:,3),'.k','markersize',15,'tag',tag) ;
    
%% LEVELSET SKELETON
    dx = norm(range(lvlst.BoundingBox,1))/50 ;
    skel = lvlst.skeleton(dx) ;
    ps = plot(skel) ;

%% POPULATE THE LEVELSET
    dx = 0.025 ;
    tag = 'population' ;
    delete(findobj(gca,'tag',tag)) ;
    P = lvlst.populate(dx,'grid',@(P)dx*(1 + abs(P(:,2))*2)) ;
    %P = [P ; lvlst.discretizeContour(dx)] ;
    P(:,end+1:3) = 0 ; % force 3D points
    c = plot3(P(:,1),P(:,2),P(:,3),'.k','markersize',8,'tag',tag) ;

%% POPULATE THE LEVELSET BOUNDARY
    dx = 0.025 ;
    tag = 'population' ;
    delete(findobj(gca,'tag',tag)) ;
    P = lvlst.populate(dx,'grid',@(P)dx*(1 + abs(P(:,2))*2),true) ;
%     P = lvlst.populate(dx,'grid',[],true) ;
    %P = [P ; lvlst.discretizeContour(dx)] ;
    P(:,end+1:3) = 0 ; % force 3D points
    c = plot3(P(:,1),P(:,2),P(:,3),'.k','markersize',8,'tag',tag) ;
    
    
%% BUILD A MESH
    delete(findobj(gca,'type','hggroup'))
    
    profile on
    mesh = pkg.geometry.mesh.distMesh(lvlst ...
                                        ...,'h0',.5 ...
                                        ... ,'pfix',lvlst.discretizeContour ...
                                        ...,'pfix',mean(lvlst.BoundingBox,1).*[.9;1.1] ...
                                        ... ,'p0',mesh.Nodes ...
                                        ... ,'debug',true ...
                                        ,'showMesh',true ...
                                        ...,'tooLongThrs',1.2 ...
                                        ...,'tooShortThrs',0.8 ...
                                        ... ,'maxCount',20 ...
                                        ...,'deltat',.75 ...
                                        ...,'qmin',0*0.005 ...
                                        ...,'Fscale',1.2 ...
                                        ... ,'bndCons',false ...
                                        ,'t_dmax',-.02 ...
                                        ,'bnd_only',lvlst.nCoord==3 ...
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
        pfix = lvlst.discretizeContour(fh) ;
        if tt==1 % Init the mesh if needed
            p0 = lvlst.populate(hmin,'iso',fh) ;
            mesh = pkg.geometry.mesh.Mesh('Nodes',p0) ;
        else
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


%% TEST RECURSIVE MESHING

dx = 30/10*lvlst.defaultDiscreteLength ;
tol = 1e-6*dx ;
maxEdgLen = 1.9*dx ;
minEdgLen = .55*dx ;
Pfix = uniquetol(lvlst.discretizeContour(dx),dx/10,'ByRows',true,'DataScale',1) ;
nPfix = size(Pfix,1) ;

% Init the mesh
Pn = [] ;
mesh = pkg.geometry.mesh.Mesh('Nodes',Pfix,'Elems',[1 2 3]) ;
clf ; axis equal ; plot(mesh)

%% Remesh
mesh.Nodes = [Pfix ; Pn] ;
mesh.Elems.NodeIdx = delaunay(mesh.Nodes) ;

% Remove faces outside the lvelset
mesh.Elems = mesh.Elems.subpart(lvlst.Function(mesh.centroid)<=-tol) ;
clf ; axis equal ; plot(mesh)

newNodes = [] ;
switch 1
    case 1 % Split interior edges
        edgToSplit = ~mesh.boundaryEdges ;
        % Do not split too short edges
        edgToSplit = edgToSplit & mesh.elemSize(mesh.Edges)>maxEdgLen ;
        newNodes = mesh.centroid(mesh.Edges.subpart(edgToSplit)) ;
    case 2 % Add a point a triangle centroid
        newNodes = mesh.centroid ;
end
% Merge too close nodes
    if ~isempty(newNodes)
    % Get nodes groups, including fixed nodes
        Pn = [Pn ; newNodes] ; 
        [~,ind] = uniquetol([Pfix ; Pn],minEdgLen,'ByRows',true,'DataScale',1,'OutputAllIndices',true) ;
    % Grouping matrix
        ii = repelem(1:numel(ind),cellfun(@numel,ind))' ;
        jj = cat(1,ind{:}) ; 
        M = sparse(ii,jj,1) ;
    % Remove nodes asociated to fixed nodes (including themselves)
        M(any(M(:,1:nPfix),2),:) = [] ;
        M(:,1:nPfix) = [] ;
    % Take the barycenter
        Pn = (M*Pn)./sum(M,2) ;
    end
% Add to mesh
clf ; axis equal ; plot(mesh)

%% To quadhex
mesh = mesh.quadhex ;
clf ; axis equal ; plot(mesh)

%% SMOOTH
mesh.LaplacianSmoothing([],100) ;
clf ; axis equal ; plot(mesh)






    



