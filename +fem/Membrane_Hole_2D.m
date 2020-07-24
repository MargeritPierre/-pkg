%% ELASTIC MEMBRANE WITH A HOLE (2D)
clc,cla%,clearvars ;

sc = 1/1 ; % units/meters
L = 0.05/sc ; % Length 
w = 0.02/sc ; % width
h = 0.03/sc ; % thickness
s = 0.00025/sc ; % approximative initial element size
Hc = [L/2 w/2] ; % center of the hole
Hr = w/2.5+eps ; % radius of the hole
El = 210e9*sc^2 ; % Young Modulus
nu = 0.3 ; % Poisson ratio
F = 1e8*s*sc ; % Applied force
Fv = [1 0] ; % Force vector
tol = s/10 ; eps*1000 ; % BC tolerance

x = linspace(0,L,ceil(L/s)+1)' ;
mesh = pkg.geometry.mesh.Mesh(x) ;
mesh.extrude([0 1]*w,ceil(w/s)) ;

lvlst = @(X)sum((X-Hc).^2,2)-Hr.^2 ;
mesh = mesh.cut(lvlst).OUT ;

%mesh.Elems = mesh.Elems.simplex ; % to tetrahedrons

mesh.LaplacianSmoothing([1 1]*0.6,10) ;

%mesh = mesh.fillHoles ;
%mesh = mesh.sortElems ;

pl = plot(mesh) ;
pl.VisibleFaces = 'none' ;
%pl.VisibleEdges = 'none' ;
%pl.HighlightBoundaryEdges = true ;
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
O = sparse(numel(ie),mesh.nNodes) ;
B = [   G{1} O ; ... E11 = dU1_dX1
        O G{2} ; ... E22 = dU2_dX2
        G{2} G{1}] ; ... 2E12 = dU1_dX2+dU2_dX1
        
% Material Stiffness Matrix
Ce = (h*El/(1-nu^2))*...
    [    1 nu 0 ; ...
         nu 1 0 ; ...
         0 0 (1-nu)/2 ; ...
    ] ;
C = repmat({sparse(Ce)},[numel(ie) 1]) ;
C = blkdiag(C{:}) ;
shiftC = (1:3:size(C,1))'+(0:2) ;
C = C(shiftC(:),shiftC(:)) ;

% Integration Weights
I = sparse(1:numel(ie),1:numel(ie),J.*W,numel(ie),numel(ie)) ;
I = repmat({I},[3 1]) ;
I = blkdiag(I{:}) ;

% Mesh Stiffness Matrix
K = B.'*I*C*B ;

% Clamped BC
BCu = abs(mesh.Nodes(:,1)-0)<=tol ;
dofBCu = repmat(BCu,[mesh.nCoord 1]) ;
K(dofBCu(:),:) = [] ;
K(:,dofBCu(:)) = [] ;
plot(mesh.Nodes(BCu,1),mesh.Nodes(BCu,2),'.r')

% Force at the other end
BCload = abs(mesh.Nodes(:,1)-L)<=tol ;
dofBCload = BCload(:) & Fv ;
f = zeros(mesh.nNodes*mesh.nCoord,1) ;
f(dofBCload(:)) = F ;
f(dofBCu) = [] ;
plot(mesh.Nodes(BCload,1),mesh.Nodes(BCload,2),'.b')

% Solve
u = zeros(mesh.nNodes,mesh.nCoord) ;
u(~dofBCu) = K\f ;

toc
profile off

% % Display
defMesh = copy(mesh) ;
% defMesh.Nodes = mesh.Nodes + u ;
% pldef = plot(defMesh) ;
% pldef.HighlightBoundaryEdges = false ;

%% 2D VISUALIZATION
EPS = mesh.interpMat(E,ie)\(reshape(B*u(:),[],3)) ;
SIG = EPS*Ce ;
devSIG = SIG - (1/2)*sum(SIG.*[1 1 0],2).*[1 1 0] ;
vmSIG = sqrt(sum((devSIG.*[1 1 2]).^2,2)/2) ;

cla
axis equal
pl = plot(mesh,'VisibleFaces','none','VisibleEdges','boundary','HighlightBoundaryEdges',false) ;
defMesh.Nodes = mesh.Nodes + u ;
pldef = plot(defMesh) ;
pldef.VisibleFaces = 'all' ;
pldef.VisibleEdges = 'none' ;
pldef.Faces.FaceColor = 'interp' ;
pldef.Faces.FaceVertexCData = EPS(:,3) ;

    
