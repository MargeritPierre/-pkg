%% COMPUTE PLANE WAVE PROPAGATION IN WAVEGUIDES

%% WAVES IN AN HOMOGENEOUS PRISMATIC WAVEGUIDE
% Parameters
    % Geometry
        L = [10 2] ; [4 .5] ; % unit cell size (mm)
        dx = min(L)/6 ; .12 ; % element size (mm)
        perVec = [] ; diag(L) ; % periodicity vectors
    % Material
        E = 70e3*(1+0.0001i) ; .5e3*(1+.00001i) ; % Young modulus (MPa)
        G = E/2.62 ; % shear modulus (MPa)
        nu = E./(2*G)-1 ;
        rho = 2700e-12 ; 1200e-12 ; % material density (tons/mm^3)
    % Wave 
        dir = [NaN NaN 1] ; % wave propagation direction [Nan: no periodicity]
        freq = linspace(100,1000e3,500) ; logspace(3,6,100) ; % wave frequency
        nModes = 50 ; % number of extracted wave modes
    % Display
        plotType = 'k' ; % 'c', 'k' , 'l' or 'g'
        logScale = false ;
        gammaMax = .1 ; inf ;
% Build the mesh
    mesh = pkg.geometry.mesh.GridMesh(L,dx) ;
    if numel(L)==2 ; mesh.setElementTypes(pkg.geometry.mesh.elements.Quad8) ; end
% Display the mesh
    clf ; axis equal tight ; view([30 30]) ;
    plot(mesh) ;
% Build the FEM matrices
    C = pkg.fem.bloch.stiffness(E,G) ;
    [K00,K0i,Kij,M,P] = pkg.fem.bloch.FEM(mesh,C,rho,[],perVec) ;
% ZERO-TH ORDER MODES
%     [U0,omega0] = pkg.fem.bloch.zeroKmodes(K00,M,nModes) ;
%     f0 = omega0/2/pi
%     plot(real(f0),'-')
% Compute the wavenumbers and modes
    [K,U,omega] = pkg.fem.bloch.solve(K00,K0i,Kij,M,freq,dir,nModes) ;
    U = reshape(P*U(:,:),[mesh.nNodes 3 nModes numel(freq)]) ;
%% DISPLAY RESULTS
    clf ;
% Theoretical wavenumbers
        b = L(1) ; h = L(2) ;
    % Longitudinal wave
        S = b*h ;
        kl = omega.*sqrt(rho/E) ;
        plot3(freq,real(kl),imag(kl),'-','linewidth',1) ;
    % Bending waves
        I11 = b^3*h/12 ;
        I22 = b*h^3/12 ;
        kb1 = sqrt(omega).*(rho*S/E/I11)^.25 ;
        kb2 = sqrt(omega).*(rho*S/E/I22)^.25 ;
        plot3(freq,real(kb1),imag(kb1),'-','linewidth',1) ;
        plot3(freq,real(kb2),imag(kb2),'-','linewidth',1) ;
    % ... with transverse shear correction
        xi = pi^2/12 ; p = rho/E ; q = rho/xi/G ;
        kb1s1 = omega.*sqrt((p+q) - sqrt((p-q)^2 + 4./(omega.^2).*(rho*S/E/I11)))/sqrt(2) ;
        plot3(freq,real(kb1s1),imag(kb1s1),'-','linewidth',1) ;
        kb1s2 = omega.*sqrt((p+q) + sqrt((p-q)^2 + 4./(omega.^2).*(rho*S/E/I11)))/sqrt(2) ;
        plot3(freq,real(kb1s2),imag(kb1s2),'-','linewidth',1) ;
        kb2s1 = omega.*sqrt((p+q) - sqrt((p-q)^2 + 4./(omega.^2).*(rho*S/E/I22)))/sqrt(2) ;
        plot3(freq,real(kb2s1),imag(kb2s1),'-','linewidth',1) ;
        kb2s2 = omega.*sqrt((p+q) + sqrt((p-q)^2 + 4./(omega.^2).*(rho*S/E/I22)))/sqrt(2) ;
        plot3(freq,real(kb2s2),imag(kb2s2),'-','linewidth',1) ;
    % Torsion wave
        a = interp1([1 1.5 2 2.5 3 4 5 6 10 1000],...
                    [.141 .196 .229 .249 .263 .281 .291 .299 .312 .333],...
                    max(b,h)/min(b,h)) ; % torsion correction, see en.wikipedia.org/wiki/Torsion_constant
        It = a*max(b,h)*min(b,h)^3 ;
        J = I11+I22 ;
        kt = omega.*(rho*J/G/It)^.5 ;
        plot3(freq,real(kt),imag(kt),'-','linewidth',1) ;
% Computed wavenumbers
    FREQ = freq + 0*K ;
    Kdir = K(:).*dir ;
    plothandle = plot3(FREQ(:),real(K(:)),imag(K(:)),'.','markersize',4,'linewidth',.1) ;
% Set plot options
    pkg.fem.bloch.setPlot(gammaMax,plotType,logScale) ;
%     pkg.fem.bloch.waveModeAnimation(mesh,Kdir,U,plothandle,true)
    

%% WAVES IN A TWO-PHASE INFINITE 3D SOLID
% Parameters
    % Geometry
        dx = .1 ; % voxel size (mm)
        L = 5*dx*[1 1 1]*(1-eps) ; % unit cell size (mm)
        subdiv = 2 ; % voxel subdivision
        perVec = diag(L) ; % periodicity vectors
    % Random phase distribution
        volfrac = .15 ;
        n = prod(ceil(L/dx)) ;
        softcell = randsample(n,floor(n*volfrac)) ; [] ; % indices
        issoft = full(sparse(softcell,1,true,n,1)) ; % logical array
    % Voxel subdivision
        issoft = reshape(issoft,ceil(L/dx)) ;
        if subdiv>1
            issoft = repelem(issoft,subdiv,subdiv,subdiv) ;
            dx = dx/subdiv ;
        end
    % Material
        E = 4e3*(1+.001i)*~issoft + 1e3*(1+.001i)*issoft ; % Young modulus (MPa)
        G = E/2.6 ; % shear modulus (MPa)
        rho = 1200e-12*~issoft + 1100e-12*issoft ; % material density (tons/mm^3)
    % Wave 
        dir = [1 0 0] ; % wave propagation direction [Nan: no periodicity]
        freq = logspace(5,7,100) ; linspace(.1e6,10e6,20) ; logspace(5,7,20) ; % wave frequency
        nModes = 10 ; % number of extracted wave modes
    % Display
        plotType = 'c' ; % 'c', 'k' , 'l' or 'g'
        logScale = true ;
        gammaMax = .1 ; inf ;
% Build the mesh
    mesh = pkg.geometry.mesh.GridMesh(L,dx) ;
% Display the mesh
    soft = copy(mesh) ; soft.Elems = mesh.Elems.subpart(issoft) ;
    clf ; axis equal tight ; view([30 30]) ;
    plot(mesh,'VisibleEdges','none','FaceAlpha',0) ;
    plot(soft,'FaceColor',[1 1 1]*.5,'FaceAlpha',0.5,'HighlightBoundaryEdges',false) ;
% Build the FEM matrices
    C = pkg.fem.bloch.stiffness(E,G) ;
    [K00,K0i,Kij,M,P] = pkg.fem.bloch.FEM(mesh,C,rho,[],perVec) ;
%% Compute the wavenumbers and modes
    [K,U,omega] = pkg.fem.bloch.solve(mesh,K00,K0i,Kij,M,P,freq,dir,nModes) ;
% DISPLAY RESULTS
    clf ;
% Theoretical wavenumbers
    rhom = mean(rho(:)) ;
    Cvoigt = mean(C,3) ;
    Creuss = 1./mean(1./C,3) ;
    % Longitudinal wave
        klr = omega.*sqrt(rhom/Creuss(1,1)) ;
        klv = omega.*sqrt(rhom/Cvoigt(1,1)) ;
        patch([freq flip(freq)],[real(klr) flip(real(klv))],[imag(klr) flip(imag(klv))],'r')
    % Transverse wave
        ktr = omega.*sqrt(rhom/Creuss(6,6)) ;
        ktv = omega.*sqrt(rhom/Cvoigt(6,6)) ;
        patch([freq flip(freq)],[real(ktr) flip(real(ktv))],[imag(ktr) flip(imag(ktv))],'b')
    set(findobj(gca,'type','patch'),'edgecolor','none','facealpha',.25)
% Computed wavenumbers
    FREQ = freq + 0*K ;
    Kdir = K(:).*dir ;
    plothandle = plot3(FREQ(:),real(K(:)),imag(K(:)),'+k','markersize',4,'linewidth',.1) ;
    %plot3(freq,real(K),imag(K),'+k','markersize',4,'linewidth',.1) ;
% Set plot options
    pkg.fem.bloch.setPlot(gammaMax,plotType,logScale) ;
    pkg.fem.bloch.waveModeAnimation(mesh,Kdir,U,plothandle) ;

    
%% WAVES IN AN HOMOGENEOUS ISOTROPIC PLATE: frequency sweep
% Parameters
    % Geometry
        L = [0 0 1] ; % unit cell size (mm)
        dx = .025 ; % element size (mm)
    % Material
        E = 210e3*(1+1i*1e-4) ; % Young modulus (MPa)
        G = E/2.6 ; % shear modulus (MPa)
        nu = E./(2*G)-1 ;
        rho = 7800e-12 ; % material density (tons/mm^3)
        sig = 0.0001*E*[0 0 0 0 0 0] ; % pre-stress (MPa) [S11 S22 S33 S13 S23 S12]
    % Wave 
        dir = [1 0 NaN] ; % wave propagation direction [Nan: no periodicity]
        freq = logspace(1,8,100) ; % wave frequency
        nModes = 80 ; % number of extracted wave modes
    % Display
        plotType = 'c' ; % 'c', 'k' , 'l' or 'g'
        logScale = true ;
        gammaMax = 1 ; inf ;
% Build the mesh
    mesh = pkg.geometry.mesh.GridMesh(L,dx) ;
% Display the mesh
    clf ; axis equal tight ; view([30 30]) ;
    plot(mesh) ;
% UNSTRESSED
    % Build the FEM matrices
        C = stiffness(E,G) ;
        [K00,K0i,Kij,M,P] = FEM(mesh,C,rho) ;
    % Compute the wavenumbers and modes
        [K0,U0] = bloch(mesh,K00,K0i,Kij,M,P,freq,dir,nModes) ;
% PRE-STRESSED
    % Build the FEM matrices
        C = pkg.fem.bloch.stiffness(E,G) ;
        [K00,K0i,Kij] = pkg.fem.bloch.FEM(mesh,C,rho,sig) ;
    % Compute the wavenumbers and modes
        [Ks,Us,omega] = pkg.fem.bloch.solve(mesh,K00,K0i,Kij,M,P,freq,dir,nModes) ;
%% DISPLAY RESULTS
    clf ;
% Theoretical wavenumbers
        h = norm(L) ;
        Q = E/(1-nu^2) ;
    % Longitudinal wave
        kl = omega.*sqrt(rho/Q) ;
        plot3(freq,real(kl),imag(kl),'-','linewidth',1) ;
    % Bending waves
        I = h^3/12 ;
        kb = sqrt(omega).*(rho*h/Q/I)^.25 ;
        plot3(freq,real(kb),imag(kb),'-','linewidth',1) ;
    % ... with transverse shear correction
        xi = pi^2/12 ; p = rho/Q ; q = rho/xi/G ;
        kbs1 = omega.*sqrt((p+q) - sqrt((p-q)^2 + 4./(omega.^2).*(rho*h/Q/I)))/sqrt(2) ;
        plot3(freq,real(kbs1),imag(kbs1),'-','linewidth',1) ;
        kbs2 = omega.*sqrt((p+q) + sqrt((p-q)^2 + 4./(omega.^2).*(rho*h/Q/I)))/sqrt(2) ;
        plot3(freq,real(kbs2),imag(kbs2),'-','linewidth',1) ;
    % Torsion wave
        kt = omega.*(rho/G)^.5 ;
        plot3(freq,real(kt),imag(kt),'-','linewidth',1) ;
% Computed wavenumbers
    FREQ = freq + 0*K0 ;
    Kdir = K0(:).*dir ;
    plothandle = plot3(FREQ(:),real(K0(:)),imag(K0(:)),'+k','markersize',4,'linewidth',.1) ;
%     plot3(freq,real(K0),imag(K0),'+k','markersize',4,'linewidth',.1) ;
%     plot3(freq,real(Ks),imag(Ks),'xr','markersize',4,'linewidth',.1) ;
% Set plot options
    pkg.fem.bloch.setPlot(gammaMax,plotType,logScale) ;
    pkg.fem.bloch.waveModeAnimation(mesh,Kdir,U0,plothandle) ;
    
    
    
%% WAVES IN AN ISOTROPIC ROUND WIRE
% Parameters
    % Geometry
        D = .3 ; % unit cell size (mm)
        dx = D/20 ; % element size (mm)
        perVec = [] ; % periodicity vectors
    % Material
        E = 210e3*(1+.00001i) ; % Young modulus (MPa)
        G = E/2.6 ; % shear modulus (MPa)
        nu = E./(2*G)-1 ;
        rho = 7800e-12 ; % material density (tons/mm^3)
    % Wave 
        dir = [NaN NaN 1] ; % wave propagation direction [Nan: no periodicity]
        freq = linspace(1e6,20e6) ; logspace(7,8,100) ; % wave frequency
        nModes = 100 ; % number of extracted wave modes
    % Display
        plotType = 'k' ; % 'c', 'k' , 'l' or 'g'
        logScale = true ;
        gammaMax = .1 ; inf ;
% Build the mesh
    mesh = pkg.geometry.levelset.Circle([0,0],D/2).mesh(dx) ;
    %if numel(L)==2 ; mesh.setElementTypes(pkg.geometry.mesh.elements.Quad8) ; end
% Display the mesh
    clf ; axis equal tight ; view([30 30]) ;
    plot(mesh) ;
% Build the FEM matrices
    C = pkg.fem.bloch.stiffness(E,G) ;
    [K00,K0i,Kij,M,P] = pkg.fem.bloch.FEM(mesh,C,rho,[],perVec) ;
% Compute the wavenumbers and modes
    [K,U,omega] = pkg.fem.bloch.solve(mesh,K00,K0i,Kij,M,P,freq,dir,nModes) ;
%% DISPLAY RESULTS
    clf ;
% Theoretical wavenumbers
    % Longitudinal wave
        S = pi*D^2/4 ;
        kl = omega.*sqrt(rho/E) ;
        plot3(freq,real(kl),imag(kl),'-','linewidth',1) ;
    % Bending waves
        I = pi*D^4/64 ;
        kb1 = sqrt(omega).*(rho*S/E/I)^.25 ;
        plot3(freq,real(kb1),imag(kb1),'-','linewidth',1) ;
    % ... with transverse shear correction
        xi = pi^2/12 ; p = rho/E ; q = rho/xi/G ;
        kb1s1 = omega.*sqrt((p+q) - sqrt((p-q)^2 + 4./(omega.^2).*(rho*S/E/I)))/sqrt(2) ;
        plot3(freq,real(kb1s1),imag(kb1s1),'-','linewidth',1) ;
        kb1s2 = omega.*sqrt((p+q) + sqrt((p-q)^2 + 4./(omega.^2).*(rho*S/E/I)))/sqrt(2) ;
        plot3(freq,real(kb1s2),imag(kb1s2),'-','linewidth',1) ;
    % Torsion wave
        kt = omega.*(rho/G)^.5 ;
        plot3(freq,real(kt),imag(kt),'-','linewidth',1) ;
% Computed wavenumbers
    FREQ = freq + 0*K ;
    Kdir = K(:).*dir ;
    plothandle = plot3(FREQ(:),real(K(:)),imag(K(:)),'+k','markersize',4,'linewidth',.1) ;
% Set plot options
    pkg.fem.bloch.setPlot(gammaMax,plotType,logScale) ;
    pkg.fem.bloch.waveModeAnimation(mesh,Kdir,U,plothandle,true)












