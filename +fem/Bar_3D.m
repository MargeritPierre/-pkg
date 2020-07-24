%% ELASTICITY OF A BAR (3D)
clc,cla%,clearvars ;

sc = 1/1 ; % units/meters
L = 0.25/sc ; % Length 
w = 0.02/sc ; % width
h = 0.03/sc ; % height
s = 0.005/sc ; % approximative initial element size
El = 210e9*sc^2 ; % Young Modulus
nu = 0.3 ; % Poisson ratio
F = 1e8*s^2*sc^2 ; % Applied force
Fv = [0 1 1] ; % force vector
tol = eps*1000 ; % BC tolerance

x = linspace(0,L,ceil(L/s)+1)' ;
mesh = pkg.geometry.mesh.Mesh(x) ;
mesh.extrude([0 1]*w,ceil(w/s)) ;
mesh.extrude([0 0 1]*h,ceil(h/s)) ;

%mesh.Elems = mesh.Elems.simplex ; % to tetrahedrons

pl = plot(mesh) ;
pl.VisibleFaces = 'none' ;
pl.VisibleEdges = 'outer' ;
pl.Edges.EdgeAlpha = 0.5 ;
pl.HighlightBoundaryEdges = false ;
%pl.ShowFrames = 'Faces' ;
axis equal
set(gca,'projection','perspective')

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
B = [   G{1} O O ; ... E11 = dU1_dX1
        O G{2} O ; ... E22 = dU2_dX2
        O O G{3} ; ... E33 = dU2_dX2
        O G{3} G{2} ; ... 2E23 = dU2_dX3+dU3_dX2
        G{3} O G{1} ; ... 2E13 = dU1_dX3+dU3_dX1
        G{2} G{1} O ] ; ... 2E12 = dU1_dX2+dU2_dX1
        
% Material Stiffness Matrix
Ce = El/(1+nu)/(1-2*nu)*...
    [    1-nu nu nu 0 0 0 ; ...
         nu 1-nu nu 0 0 0 ; ...
         nu nu 1-nu 0 0 0 ; ...
         0 0 0 (1-2*nu)/2 0 0 ; ...
         0 0 0 0 (1-2*nu)/2 0 ; ...
         0 0 0 0 0 (1-2*nu)/2 ; ...
    ] ;
C = repmat({sparse(Ce)},[numel(ie) 1]) ;
C = blkdiag(C{:}) ;
shiftC = (1:6:size(C,1))'+(0:5) ;
C = C(shiftC(:),shiftC(:)) ;

% Integration Weights
I = sparse(1:numel(ie),1:numel(ie),J.*W,numel(ie),numel(ie)) ;
I = repmat({I},[6 1]) ;
I = blkdiag(I{:}) ;

% Mesh Stiffness Matrix
K = B.'*I*C*B ;

% Clamped BC
BCu = abs(mesh.Nodes(:,1)-0)<=tol ;
dofBCu = repmat(BCu,[mesh.nCoord 1]) ;
K(dofBCu(:),:) = [] ;
K(:,dofBCu(:)) = [] ;
plot3(mesh.Nodes(BCu,1),mesh.Nodes(BCu,2),mesh.Nodes(BCu,3),'.r')

% Force at the other end in the Z direction
BCload = abs(mesh.Nodes(:,1)-L)<=tol ;
dofBCload = BCload(:) & Fv ;
f = zeros(mesh.nNodes*mesh.nCoord,1) ;
f(dofBCload(:)) = F ;
f(dofBCu) = [] ;
plot3(mesh.Nodes(BCload,1),mesh.Nodes(BCload,2),mesh.Nodes(BCload,3),'.b')

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


%% 3D VISUALIZATION
EPS = mesh.interpMat(E,ie)\(reshape(B*u(:),[],6)) ;
SIG = EPS*Ce ;
devSIG = SIG - (1/3)*sum(SIG.*[1 1 1 0 0 0],2).*[1 1 1 0 0 0] ;
vmSIG = sqrt(sum((devSIG.*[1 1 1 6 6 6]).^2,2)/2) ;

cla
axis equal
pl = plot(mesh,'VisibleFaces','none','VisibleEdges','boundary','HighlightBoundaryEdges',false) ;
defMesh.Nodes = mesh.Nodes + u ;
pldef = plot(defMesh) ;
pldef.VisibleFaces = 'all' ;
pldef.VisibleEdges = 'none' ;
pldef.Faces.FaceColor = 'interp' ;
pldef.Faces.FaceVertexCData = vmSIG ;