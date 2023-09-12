%% WAVES IN AN HOMOGENEOUS PRISMATIC WAVEGUIDE
% Parameters
    % Geometry
        L = [4 .5] ; % unit cell size (mm)
        dx = .12 ; % element size (mm)
        perVec = [] ; % no periodicity conditions
    % Material
        E = .5e3*(1+.00000i) ; % Young modulus (MPa)
        G = E/2.9 ; % shear modulus (MPa)
        nu = E./(2*G)-1 ;
        rho = 1200e-12 ; % material density (tons/mm^3)
    % Wave 
        dir = [NaN NaN 1] ; % wave propagation direction [Nan: no periodicity]
        freq = logspace(3,6,100) ; % wave frequency
        nModes = 50 ; % number of extracted wave modes
    % Display
        plotType = 'k' ; % 'c', 'k' , 'l' or 'g'
        logScale = true ;
        gammaMax = .1 ; inf ;
% Build the mesh
    mesh = pkg.geometry.mesh.GridMesh(L,dx) ;
    %if numel(L)==2 ; mesh.setElementTypes(pkg.geometry.mesh.elements.Quad8) ; end
% Display the mesh
    clf ; axis equal tight ; view([30 30]) ;
    plot(mesh) ;
% Build the FEM matrices
    C = pkg.fem.bloch.stiffness(E,G) ;
    [K00,K0i,Kij,M,P] = pkg.fem.bloch.FEM(mesh,C,rho,[],[]) ;
% Compute the wavenumbers and modes
    [K,U,omega] = pkg.fem.bloch.solve(K00,K0i,Kij,M,freq,dir,nModes) ;
    U = reshape(P*U(:,:),[mesh.nNodes 3 size(U,2:ndims(U))]) ;
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
    plothandle = plot3(FREQ(:),real(K(:)),imag(K(:)),'+k','markersize',4,'linewidth',.1) ;
% Set plot options
    pkg.fem.bloch.setPlot(gammaMax,plotType,logScale) ;
    pkg.fem.bloch.waveModeAnimation(mesh,Kdir,U,plothandle,false)
    
%% SOLVE FOR ZERO MODES
nModes0 = 50 ;
    [U0,omega0] = pkg.fem.bloch.zeroKmodes(K00,M,nModes0) ;
    U0 = reshape(P*U0(:,:),[mesh.nNodes 3 size(U0,2:ndims(U0))]) ;
    f0 = omega0/2/pi 

%% PLOT ZERO MODES
ff = 6 ;
    uu = U0(:,:,ff) ;
    uu = .025*uu./max(abs(uu(:))).*max(range(mesh.Nodes,1)) ;
    clf ; axis equal
    pl0 = plot(mesh,'FaceColor','none') ;
    pl = plot(mesh,'Deformation',uu) ;
    view([50 30])
    title("$f^c_{"+string(ff)+"}$ = "+string(f0(ff))+"Hz")

%% SOLVE WITH ZERO MODES AS A BASIS
% Project on the basis
    Pz = reshape(U0(:,:,1:20),[],20) ;  Ua(:,1:20) ;
    Mz = Pz'*M*Pz ;
    K00z = Pz'*K00*Pz ;
    for ii = 1:numel(K0i) ; K0iz{ii} = Pz'*K0i{ii}*Pz ; end
    for ii = 1:numel(Kij) ; Kijz{ii} = Pz'*Kij{ii}*Pz ; end
% Solve
    [Kz,Uz] = pkg.fem.bloch.solve(K00z,K0iz,Kijz,Mz,freq,dir) ;
    Uz = reshape(Pz*Uz(:,:),[mesh.nNodes 3 size(Uz,2:ndims(Uz))]) ;
% Display
    clf ;
% Reference wavenumbers
    FREQ = freq + 0*K ;
    Kdir = K(:).*dir ;
    plothandle = plot3(FREQ(:),real(K(:)),imag(K(:)),'+k','markersize',4,'linewidth',.1) ;
% Approximation
    FREQz = freq + 0*Kz ;
    Kdirz = Kz(:).*dir ;
    plotzhandle = plot3(FREQz(:),real(Kz(:)),imag(Kz(:)),'or','markersize',4,'linewidth',.1) ;
% Set plot options
    pkg.fem.bloch.setPlot(gammaMax,plotType,logScale) ;
%     pkg.fem.bloch.waveModeAnimation(mesh,Kdir,U,plothandle,false)

%% MODEL REDUCTION VIA POD
zeroMean = true ;
uu = U ;
if zeroMean
    uu = reshape(uu,mesh.nNodes,3,[]) ;
    uu = uu-mean(uu,1) ;
end
uu = reshape(uu,mesh.nNodes*3,[]) ;
% uu = uu.*abs(K(:)') ;
FREQ = freq + 0*K ;
uu = uu./abs(FREQ(:)') ;
uu(:,abs(imag(K)./real(K))>gammaMax) = [] ;
[Ua,siga] = eig(uu*uu','vector') ;
[siga,isort] = sort(siga,'descend') ;
Ua = Ua(:,isort) ;
if zeroMean
    Ua = [kron(eye(3),ones(mesh.nNodes,1)) Ua] ;
end
clf ; plot(1-sqrt(cumsum(siga.^2)/sum(siga.^2))) ; set(gca,'yscale','log') ;

%% PLOT POD MODES
mm = 8 ;
    uu = reshape(Ua(:,mm),mesh.nNodes,3) ;
    uu = real(uu) ; 
    uu = .025*uu./max(abs(uu(:))).*max(range(mesh.Nodes,1)) ;
    clf ; axis equal
    pl0 = plot(mesh,'FaceColor','none') ;
    pl = plot(mesh,'Deformation',uu) ;
    view([50 30])





    