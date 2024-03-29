clc
clearvars

%% RANDOM P1-SIMPLEX MESH GENERATION

clf
    axis equal
    axis ij
    set(gca,'xtick',[],'ytick',[]) ;
    axis off

tol = 0.02 ;
Nodes = rand(1000000,2) ;
Nodes = uniquetol(Nodes,tol,'ByRows',true) ;
Tris = delaunay(Nodes(:,1:2)) ;

tic
mesh = pkg.geometry.mesh.Mesh('Nodes',Nodes,'Elems',Tris) ;
h = plot(mesh) ;
toc


%% QUAD MESH ON A GRID
N = 10 ;
[X,Y] = meshgrid(linspace(0,1,N),linspace(0,1,N)) ;
Elems = [1:N-1 ; 2:N ; N+2:2*N ; N+1:2*N-1]' ;
Elems = repmat(Elems,[N-1 1]) + N*kron((0:N-2)',Elems*0+1) ;
Nodes = [X(:) Y(:)] ;
Elems = Elems(:,[1 2 4 3]) ; % for lagrange quad compatibility

mesh = pkg.geometry.mesh.Mesh('Nodes',Nodes,'Elems',Elems) ;

cla
plot(mesh)
axis equal


%% LAPLACIAN SMOOTHING

lmbda = [5/10 0/10] ;
iterations = 1000 ;
tol = 1e-9 ;
mesh.sortElems ;

tic
mesh2 = mesh.LaplacianSmoothing(lmbda,iterations,tol) ;
toc
cla ;
h = plot(mesh2) ;


%% CATMULL-CLARK SUBDIVISION

iterations = 1 ;

tic
mesh = mesh.CatmullClark(iterations) ;
toc
delete(h) ;
h = plot(mesh) ;


%% ELEMENT ORDER CHANGE

[mesh.nNodes mesh.nElems]
mesh.setElementOrder(4) ;
[mesh.nNodes mesh.nElems]

clf
h2 = plot(mesh) ;
% h2.Patches.EdgeColor = 'r' ;
% h2.Lines.Color = 'r' ;
axis equal


%% MESH TRANSFORMATION & DUPLICATION
clc

% N = 10 ;
% [X,Y] = meshgrid(linspace(0,1,N),linspace(0,1,N)) ;
% Elems = [1:N-1 ; 2:N ; N+2:2*N ; N+1:2*N-1]' ;
% Elems = repmat(Elems,[N-1 1]) + N*kron((0:N-2)',Elems*0+1) ;
% Nodes = [X(:) Y(:)] ;
% Elems = Elems(:,[1 2 4 3]) ; % for lagrange quad compatibility
% 
% mesh = pkg.geometry.mesh.Mesh('Nodes',Nodes,'Elems',Elems) ;

clf ;
h = plot(mesh) ;
h.Patches.EdgeColor = 'r' ;
h.Lines.Color = 'r' ;
axis equal

mesh.rotate([10 20 30 45]*pi/180,[1 1]*0) ;

h2 = plot(mesh) ;

%% POINTS INSIDE AN ELEMENT

elmt = pkg.geometry.mesh.elements.LagrangeElement('tet',2) ;

bbox = elmt.localCoordinatesDomain ;
bbox = bbox + 0.5*range(bbox,1).*[-1; 1] ; % extend bbox
E = (rand(1000000,elmt.nDims).*range(bbox,1))+bbox(1,:) ;
tol = 1e-9 ;

tic ;
[in,on] = elmt.isInside(E,tol) ;
toc

clf ;
h = plot(elmt) ;
pIn = patch('vertices',E,'faces',find(in & ~on),'facecolor','none','edgecolor','none','marker','.','markeredgecolor','r','markersize',5) ;
pOn = patch('vertices',E,'faces',find(on),'facecolor','none','edgecolor','none','marker','.','markeredgecolor','b','markersize',5) ;
%pOut = patch('vertices',E,'faces',find(~in),'facecolor','none','edgecolor','none','marker','.','markeredgecolor','g','markersize',5) ;
axis equal


%% POINTS INSIDE THE MESH
clc

pts = rand(100,mesh.nCoord) ;

in = mesh.isInside(pts) ; % Which point is on the mesh ?
M = mesh.isInside(pts,mesh.Elems) ; % Which point is on which mesh element ?

if exist('p','var') ; delete(p) ; end
p = gobjects(0) ;
if any(in) ; p(end+1) = plot(pts(in,1),pts(in,2),'.r','markersize',10) ; end
if any(~in) ; p(end+1) = plot(pts(~in,1),pts(~in,2),'.b','markersize',10) ; end

%% POINT LOCALIZATION

mesh.setElementOrder(2) ;
P = rand(10000,mesh.nCoord) ; 

clf ;
h = plot(mesh) ; 
axis equal

%E = mesh.localize(P,mesh.Elems,mesh.X.Values,mesh.defaultTolerance) ;

%profile on
tic ; 
M = mesh.isInside(P,mesh.Elems) ; 
[pp,ee] = find(M) ;  
E = mesh.localize(P,mesh.Elems,mesh.X.Values,mesh.defaultTolerance,pp,ee) ;
toc
%profile viewer

    
%% IMPORT FROM navDIC

global hd
Seed = hd.Seeds(1) ;
mesh = pkg.geometry.mesh.Mesh('Nodes',Seed.Points,'Elems',Seed.Triangles) ;

clf 
    h = mesh.plot() ;
    axis equal
    axis off
    axis ij
    
    





