%% WAVES IN A TWO-PHASE MATRIX-INCLUSION 2D UNIT CELL
clc,clearvars
% Parameters
    % Geometry
        L = .6770*[1 1] ; % unit cell size (mm)
        a = .3385*[1 1] ; % square inclusion size (mm)
        dx = min(a)/8 ; % approximate element edge length (mm)
        elemType = 'quad' ;
        subdiv = 0 ;
        perVec = diag(L) ; % periodicity vectors
    % Material properties
        C11_hard = 7.1e3 ; C44_hard = 1.5e3 ; % Moduli (MPa)
        rho_hard = 1177e-12 ; % material density (tons/mm^3)
        C11_soft = 3.8e3 ; C44_soft = 0.5e3 ; % Moduli (MPa)
        rho_soft = 1104e-12 ; % material density (tons/mm^3)
        %C11_soft = C11_hard ; C44_soft = C44_hard ; rho_soft = rho_hard ; 
    % Wave 
        dir = [1 0 0] ; % wave propagation direction [Nan: no periodicity]
        freq = linspace(10000,3e6,300) ; % wave frequency (Hz)
        nModes = 20 ; % number of extracted wave modes
    % Display
        plotType = 'k' ; % 'c', 'k' , 'l' or 'g'
        logScale = false ;
        gammaMax = .1 ; inf ;
% Meshing
    soft = pkg.geometry.levelset.Rectangle(.5*a.*[-1;1]) ; % soft phase
    switch elemType
        case 'tri'
            hard = pkg.geometry.levelset.Rectangle(.5*L.*[-1;1]) ; % hard phase
            hard = (hard - soft) ; % boolean operation
            mesh = merge(hard.mesh(dx),soft.mesh(dx)) ; % mesh the two phases independently then merge the meshes
        otherwise % 'quad' by default
            mesh =  pkg.geometry.mesh.GridMesh(a/2,dx) ... % central part
                    + pkg.geometry.mesh.GridMesh((L-a)/2,dx).move(a/2) ... top-right corner
                    + pkg.geometry.mesh.GridMesh([a(1) L(2)-a(2)]/2,dx).move([0 1].*a/2) ... top row
                    + pkg.geometry.mesh.GridMesh([L(1)-a(1) a(2)]/2,dx).move([1 0].*a/2) ... right column
                    ;
            %mesh = mesh.rotate(pi/2*(0:3)) ;
            mesh = mesh.replicate(mesh.Nodes.*cat(3,[1 1],[-1 1],[-1 -1],[1 -1]),false) ;
            mesh.clean() ;
            mesh.sortElems ;
    end
    mesh.CatmullClark(subdiv) ;
    if strcmp(elemType,'Quad8')
        mesh.setElementTypes(pkg.geometry.mesh.elements.Quad8) ;
    end
% Local material assignment for each element
    issoft = soft.inside(mesh.centroid) ; % is the element centroid inside the soft phase ?
    G = C44_soft.*issoft + C44_hard.*~issoft; % shear modulus (MPa)
    lmbda = (C11_soft.*issoft + C11_hard.*~issoft) - 2*G ; % Lamé coeff. lambda (MPa)
    E = G.*(3.*lmbda+2.*G)./(lmbda+G) ; % Young modulus (MPa)
    rho = rho_soft.*issoft + rho_hard.*~issoft ; % material density (tons/mm^3)
% Mesh display
    clf ; axis equal tight ; plot(mesh,'CData',[1 1 1].*(~issoft)+.65,'VisibleNodes','all') ;
% Build the FEM matrices
    C = pkg.fem.bloch.stiffness(E,G) ; % material stiffnesses [6 6 mesh.nElems] 
    [K00,K0i,Kij,M,P] = pkg.fem.bloch.FEM(mesh,C,rho,[],perVec) ;
%% Compute the wavenumbers and modes
    [K,U,omega] = pkg.fem.bloch.solve(K00,K0i,Kij,M,freq,dir,nModes) ;
    U = reshape(P*U(:,:),[mesh.nNodes 3 size(U,2:ndims(U))]) ;
%% DISPLAY RESULTS
    clf ;
        plotType = 'k' ; % 'c', 'k' , 'l' or 'g'
        logScale = false ;
        gammaMax = .1 ; inf ;
% Theoretical wavenumbers
    if 0
        cl_hard = sqrt(C11_hard./rho_hard) ;
        ct_hard = sqrt(C44_hard./rho_hard) ;
        cl_soft = sqrt(C11_soft./rho_soft) ;
        ct_soft = sqrt(C44_soft./rho_soft) ;
        plot3(freq,real(omega./cl_hard),imag(omega./cl_hard),'-','linewidth',1) ;
        plot3(freq,real(omega./ct_hard),imag(omega./ct_hard),'-','linewidth',1) ;
        plot3(freq,real(omega./cl_soft),imag(omega./cl_soft),'-','linewidth',1) ;
        plot3(freq,real(omega./ct_soft),imag(omega./ct_soft),'-','linewidth',1) ;
    end
% Computed wavenumbers
    FREQ = freq + 0*K ;
    Kdir = K(:).*dir ;
    plothandle = plot3(FREQ(:),real(K(:)),imag(K(:)),'.k','markersize',4,'linewidth',.1) ;
    %plot3(freq,real(K),imag(K),'+k','markersize',4,'linewidth',.1) ;
% Set plot options
    pkg.fem.bloch.setPlot(gammaMax,plotType,logScale) ;
    %pkg.fem.bloch.waveModeAnimation(mesh,Kdir,U,plothandle,false) ;