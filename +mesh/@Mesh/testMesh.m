clc
clearvars

%% RANDOM MESH GENERATION

clf
    axis equal
    axis ij
    set(gca,'xtick',[],'ytick',[]) ;
    axis off

tol = 0.05 ;
Nodes = rand(10000,2) ;
Nodes = uniquetol(Nodes,tol,'ByRows',true) ;
Tris = delaunay(Nodes(:,1:2)) ;

tic
mesh = pkg.mesh.Mesh('Nodes',Nodes,'Elems',Tris) ;
h = plot(mesh) ;
toc


%% LAPLACIAN SMOOTHING

lmbda = [1 0.01] ;
iterations = 1000 ;
tol = 1e-9 ;

tic
mesh.LaplacianSmoothing(lmbda,iterations,tol) ;
toc
delete(h) ;
h = plot(mesh) ;


%% CATMULL-CLARK SUBDIVISION

iterations = 1 ;

tic
mesh = mesh.CatmullClark(iterations) ;
toc
delete(h) ;
h = plot(mesh) ;


%% MESH TRANSFORMATION
clc

N = 10 ;
[X,Y] = meshgrid(linspace(0,1,N),linspace(0,1,N)) ;
Elems = [1:N-1 ; 2:N ; N+2:2*N ; N+1:2*N-1]' ;
Elems = repmat(Elems,[N-1 1]) + N*kron((0:N-2)',Elems*0+1) ;
Nodes = [X(:) Y(:)] ;
Elems = Elems(:,[1 2 4 3]) ; % for lagrange quad compatibility

mesh = pkg.mesh.Mesh('Nodes',Nodes,'Elems',Elems) ;

clf ;
h = plot(mesh) ;
h.Patches.EdgeColor = 'r' ;
h.Lines.Color = 'r' ;
axis equal

mesh.rotate([10 20 30 45]*pi/180,[1 1]*0) ;

h2 = plot(mesh) ;

    
%% IMPORT FROM navDIC

global hd
Seed = hd.Seeds(1) ;
mesh = pkg.mesh.Mesh('Nodes',Seed.Points,'Elems',Seed.Triangles) ;

clf 
    h = mesh.plot() ;
    axis equal
    axis off
    axis ij
    
    





