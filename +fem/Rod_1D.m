%% ELASTICITY OF A ROD (1D)
clc,cla%,clearvars ;

sc = 1 ; % units/meters
L = 0.5/sc ; % Length 
w = 0.01/sc ; % width
h = 0.04/sc ; % height
s = 0.02/sc ; % approximative initial element size
El = 210e9*sc^2 ; % Young Modulus
nu = 0.3 ; % Poisson ratio
F = 10000000000*s ; % Applied force
Fv = [1] ; % Force vector
tol = eps*1000 ; % BC tolerance

x = linspace(0,L,ceil(L/s)+1)' ;
mesh = pkg.geometry.mesh.Mesh(x) ;

pl = plot(mesh) ;
%pl.VisibleFaces = 'none' ;
%pl.VisibleEdges = 'none' ;
%pl.HighlightBoundaryEdges = false ;
%pl.ShowFrames = 'Faces' ;
axis equal
set(gca,'projection','orthographic')

profile on
tic

% Quadrature locations
[E,iep] = mesh.Elems.getListOf('GaussIntegrationPoints') ;
[W,iew] = mesh.Elems.getListOf('GaussIntegrationWeights') ;
if ~isequal(iep,iew) ; error('Quadrature points and weights do not match') ; end
ie = iep ;

% Quadrature in the reference element
J = mesh.detJacobian(E,ie) ;

% Gradient matrix
G = mesh.gradMat(E,ie) ;

% Strain computation
B = G{1} ; % E11 = dU1_dX1

% Integration Weights
I = sparse(1:numel(ie),1:numel(ie),J.*W,numel(ie),numel(ie)) ;

% Mesh Stiffness Matrix
K = (El*(w*h))*B.'*I*B ;

% Clamped BC
BCu = abs(mesh.Nodes(:,1)-0)<=tol ;
dofBCu = repmat(BCu,[mesh.nCoord 1]) ;
K(dofBCu(:),:) = [] ;
K(:,dofBCu(:)) = [] ;
plot(mesh.Nodes(BCu,1),mesh.Nodes(BCu,1)*0,'.r')

% Force at the other end
BCload = abs(mesh.Nodes(:,1)-L)<=tol ;
dofBCload = BCload(:) & Fv ;
f = zeros(mesh.nNodes*mesh.nCoord,1) ;
f(dofBCload(:)) = F ;
f(dofBCu) = [] ;
plot(mesh.Nodes(BCload,1),mesh.Nodes(BCload,1)*0,'.b')

% Solve
u = zeros(mesh.nNodes,mesh.nCoord) ;
u(~dofBCu) = K\f ;

toc
profile off

% Display
defMesh = copy(mesh) ;
defMesh.Nodes = mesh.Nodes + u ;
pldef = plot(defMesh) ;
pldef.HighlightBoundaryEdges = false ;