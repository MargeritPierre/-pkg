%% BUCKLING BEAM: ITERATIVE SCHEME

%% EULER BUCKLING
clc,cla%,clearvars ;

sc = 1/1 ; % units/meters
L = 1/sc ; % beam length 
w = 0.02/sc ; % beam width
h = 0.03/sc ; % beam thickness
s = 0.0025/sc ; % approximative initial element size
El = 210e9*sc^2 ; % Young Modulus
nu = 0.3 ; % Poisson ratio
F1 = -1e6*s*sc ; % Incrementally applied force
F2 = F1*1e-3 ; % Force perturbation in transverse direction
Nit = 200 ; % number of iterations
tol = s/10 ; eps*1000 ; % BC tolerance

x = linspace(0,L,ceil(L/s)+1)' ;
mesh = pkg.geometry.mesh.Mesh(x) ;
mesh.move([0 -1]*w/2) ;
mesh.extrude([0 1]*w,ceil(w/s)) ;

% Init Display
pl = plot(mesh) ;
pl.VisibleFaces = 'all' ;
pl.Faces.FaceColor = 'interp' ;
%pl.VisibleEdges = 'none' ;
%pl.HighlightBoundaryEdges = true ;
%pl.ShowFrames = 'Faces' ;
axis equal
set(gca,'projection','orthographic')

% CONSTANT OBJECTS

    % Quadrature locations (in elmt local coordinates)
    [E,iep] = mesh.Elems.getListOf('GaussIntegrationPoints') ;
    [W,iew] = mesh.Elems.getListOf('GaussIntegrationWeights') ;
    if ~isequal(iep,iew) ; error('Quadrature points and weights do not match') ; end
    ie = iep ;

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

    % Quadrature in the reference element
    J = mesh.detJacobian(E,ie) ;
    
    % Boundary Conditions
    % Clamped end
        BCu = abs(mesh.Nodes(:,1)-0)<=tol ;
        dofBCu = repmat(BCu,[1 mesh.nCoord]) ;
    % Applied load at othe end (non-following)
        BCload = abs(mesh.Nodes(:,1)-L)<=tol ;
        f = zeros(mesh.nNodes,mesh.nCoord) ;
        f(BCload(:),1) = F1 ;
        f(BCload(:),2) = F2 ;
        f(dofBCu(:)) = [] ;
    
    
% ITERATIVE SCHEME
    SIG = zeros(numel(ie),3) ;
    step = linspace(0,1,Nit) ;
    for ii = 1:Nit
        t = tic ;
        
        % Gradient matrix
        G = mesh.gradMat(E,ie) ;

        % Strain computation
        O = sparse(numel(ie),mesh.nNodes) ;
        B = [   G{1} O ; ... E11 = dU1_dX1
                O G{2} ; ... E22 = dU2_dX2
                G{2} G{1}] ; ... 2E12 = dU1_dX2+dU2_dX1
        

        % Integration Weights
        I = sparse(1:numel(ie),1:numel(ie),J.*W,numel(ie),numel(ie)) ;
        I = repmat({I},[3 1]) ;
        I = blkdiag(I{:}) ;

        % Mesh Stiffness Matrix
        K = B'*I*C*B ;

        % Clamped BC
        K(dofBCu(:),:) = [] ;
        K(:,dofBCu(:)) = [] ;
        
        % Solve
        u = zeros(mesh.nNodes,mesh.nCoord) ;
        u(~dofBCu(:)) = K\(f(:)*step(ii)-B(:,~dofBCu(:))'*I*SIG(:)) ;
        
        % Compute stresses
        SIG = SIG + reshape(C*B*u(:),[],3) ;
        
        toc(t)
        
        % Display
        mesh.Nodes = mesh.Nodes + u ;
        pl.update ;
        pl.Faces.FaceVertexCData = sqrt(sum(u.^2,2)) ;
        drawnow ;
    end


%% SNAP-THROUGH
% clc,cla%,clearvars ;
% 
% sc = 1/1 ; % units/meters
% L = 1/sc ; % beam length 
% w = 0.01/sc ; % beam width
% h = 0.03/sc ; % beam thickness
% a = L/10 ; % shape amplitude
% s = 0.0025/sc ; % approximative initial element size
% El = 210e9*sc^2 ; % Young Modulus
% nu = 0.3 ; % Poisson ratio
% F2 = -2e7*s*sc ; % Max. force in transverse direction
% F1 = 0*F2*1e-3 ; % Max. force in longitudinal direction
% Nit = 200 ; % number of iterations
% tol = s/10 ; % BC tolerance
% 
% x = linspace(0,L,ceil(L/s)+1)' ;
% x = [x(:) a*(1-cos(2*pi*x(:)/L))] ;
% mesh = pkg.geometry.mesh.Mesh(x) ;
% 
% mesh.offset(linspace(-1,1,ceil(w/s)+1)*w/2,true) ;
% 
% % Init Display
% pl = plot(mesh) ;
% pl.VisibleFaces = 'all' ;
% pl.Faces.FaceColor = 'interp' ;
% %pl.VisibleEdges = 'none' ;
% %pl.HighlightBoundaryEdges = true ;
% %pl.ShowFrames = 'Faces' ;
% axis equal
% set(gca,'projection','orthographic')
% 
% % CONSTANT OBJECTS
% 
%     % Quadrature locations (in elmt local coordinates)
%     [E,iep] = mesh.Elems.getListOf('GaussIntegrationPoints') ;
%     [W,iew] = mesh.Elems.getListOf('GaussIntegrationWeights') ;
%     if ~isequal(iep,iew) ; error('Quadrature points and weights do not match') ; end
%     ie = iep ;
% 
%     % Material Stiffness Matrix
%     Ce = (h*El/(1-nu^2))*...
%         [    1 nu 0 ; ...
%              nu 1 0 ; ...
%              0 0 (1-nu)/2 ; ...
%         ] ;
%     C = repmat({sparse(Ce)},[numel(ie) 1]) ;
%     C = blkdiag(C{:}) ;
%     shiftC = (1:3:size(C,1))'+(0:2) ;
%     C = C(shiftC(:),shiftC(:)) ;
% 
%     % Quadrature in the reference element
%     J = mesh.detJacobian(E,ie) ;
%     
%     % Boundary Conditions
%     % Clamped end
%         BCu = abs(mesh.Nodes(:,1)-0)<=tol | abs(mesh.Nodes(:,1)-L)<=tol ;
%         dofBCu = repmat(BCu,[1 mesh.nCoord]) ;
%     % Applied load at othe end (non-following)
%         BCload = abs(mesh.Nodes(:,1)-L/2)<=tol ;
%         f = zeros(mesh.nNodes,mesh.nCoord) ;
%         f(BCload(:),1) = F1 ;
%         f(BCload(:),2) = F2 ;
%         f(dofBCu(:)) = [] ;
%     
%     
% % ITERATIVE SCHEME
%     SIG = zeros(numel(ie),3) ;
%     step = linspace(0,1,Nit) ;
%     for ii = 1:Nit
%         t = tic ;
%         
%         % Gradient matrix
%         G = mesh.gradMat(E,ie) ;
% 
%         % Strain computation
%         O = sparse(numel(ie),mesh.nNodes) ;
%         B = [   G{1} O ; ... E11 = dU1_dX1
%                 O G{2} ; ... E22 = dU2_dX2
%                 G{2} G{1}] ; ... 2E12 = dU1_dX2+dU2_dX1
%         
% 
%         % Integration Weights
%         I = sparse(1:numel(ie),1:numel(ie),J.*W,numel(ie),numel(ie)) ;
%         I = repmat({I},[3 1]) ;
%         I = blkdiag(I{:}) ;
% 
%         % Mesh Stiffness Matrix
%         K = B'*I*C*B ;
% 
%         % Clamped BC
%         K(dofBCu(:),:) = [] ;
%         K(:,dofBCu(:)) = [] ;
%         
%         % Solve
%         u = zeros(mesh.nNodes,mesh.nCoord) ;
%         u(~dofBCu(:)) = K\(f(:)*step(ii)-B(:,~dofBCu(:))'*I*SIG(:)) ;
%         
%         % Compute stresses
%         SIG = SIG + reshape(C*B*u(:),[],3) ;
%         
%         toc(t)
%         
%         % Display
%         mesh.Nodes = mesh.Nodes + u ;
%         pl.update ;
%         pl.Faces.FaceVertexCData = sqrt(sum(u.^2,2)) ;
%         drawnow ;
%     end
