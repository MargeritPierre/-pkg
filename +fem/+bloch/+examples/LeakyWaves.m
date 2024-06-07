%% WAVES IN AN HOMOGENEOUS ISOTROPIC PLATE: frequency sweep
clc,clear all
% Parameters
    planeProblem = true ;
    % Plate
        hp = 10 ; % plate thickness
        Ep = 70e3*(1+0i*1e-4) ; % Young modulus (MPa)
        Gp = Ep/2.66 ; % shear modulus (MPa)
        rhoP = 2600e-12 ; % material density (tons/mm^3)
        nuP = Ep./(2*Gp)-1 ;
    % Fluid
        hf = 2 ;
        cf = 1500e3 ; % wave velocity in the fluid (mm/s)
        %omega0 = 2*pi*100e3 ; % "reference frequency (rad/s)
        %muF = .7e-6 ; % fluid viscosity at ref freq (MPa.s)
        rhoF = 1200e-12 ; % fluid density (tons/mm^3)
        Gf = 3500 ; %-1i*muF*omega0 ; % shear "Moduli" of the fluid (MPa)
        Mf = cf.^2*rhoF ;
        Ef = Gf*(3*Mf-4*Gf)/(Mf-Gf) ; % see https://fr.wikipedia.org/wiki/Coefficient_de_Lamé
        %Ef = 3.5e3 ; Gf = Ef/2.66 ;
        %rhoF = rhoP ; Ef = Ep ; Gf = Gp ;
    % Mesh discretization
        dx = min([hp,hf])/30 ; % element size (mm)
    % Wave 
        dir = [1 0 NaN] ; % wave propagation direction [Nan: no periodicity]
        freq = linspace(1e-0,1e6,100) ; % wave frequency
        nModes = 50 ; % number of extracted wave modes
    % Display
        plotType = 'k' ; % 'c', 'k' , 'l' or 'g'
        logScale = false ;
        gammaMax = 1 ; inf ;
% Build the mesh
    mesh = pkg.geometry.mesh.GridMesh(hp*[0 1]'.*[0 0 1],dx) ;
    mesh = mesh.merge(pkg.geometry.mesh.GridMesh((hp + hf.*[0 1]').*[0 0 1],dx)) ;
% Display the mesh
    clf ; axis equal tight ; view([30 30]) ;
    plot(mesh) ;
% Material properties
    [ee,we,ie] = mesh.integration() ;
    ze = mesh.interpMat(ee,ie)*mesh.Nodes(:,3) ;
    isInPlate = ze<=hp ;
    E = Ep.*isInPlate + Ef.*~isInPlate ;
    G = Gp.*isInPlate + Gf.*~isInPlate ;
    rho = rhoP.*isInPlate + rhoF.*~isInPlate ;
% Build the FEM matrices
    C = pkg.fem.bloch.stiffness(E,G) ;
    [K00,K0i,Kij,M,P] = pkg.fem.bloch.FEM(mesh,C,rho,[],[]) ;
    if planeProblem % remove the u2 component from the problem
        dofs = true(mesh.nNodes,1) & [1 0 1] ; % take only the components u1 and u3
        P = P(dofs,dofs) ; 
        M = M(dofs,dofs) ;
        K00 = K00(dofs,dofs) ;
        K0i = cellfun(@(K)K(dofs,dofs),K0i,'uni',false) ;
        Kij = cellfun(@(K)K(dofs,dofs),Kij,'uni',false) ;
    end
% Compute the wavenumbers and modes
    [K0,U0,omega] = pkg.fem.bloch.solve(K00,K0i,Kij,M,freq,dir,nModes) ;
% DISPLAY RESULTS
    clf ;
% Theoretical wavenumbers
        Qp = Ep/(1-nuP^2) ;
    % Longitudinal wave
        kl = omega.*sqrt(rhoP/Qp) ;
        plot3(freq,real(kl),imag(kl),'-','linewidth',1) ;
    % Bending waves
        I = hp^3/12 ;
        kb = sqrt(omega).*(rhoP*hp/Qp/I)^.25 ;
        plot3(freq,real(kb),imag(kb),'-','linewidth',1) ;
    % ... with transverse shear correction
        xi = pi^2/12 ; p = rhoP/Qp ; q = rhoP/xi/Gp ;
        kbs1 = omega.*sqrt((p+q) - sqrt((p-q)^2 + 4./(omega.^2).*(rhoP*hp/Qp/I)))/sqrt(2) ;
        plot3(freq,real(kbs1),imag(kbs1),'-','linewidth',1) ;
        kbs2 = omega.*sqrt((p+q) + sqrt((p-q)^2 + 4./(omega.^2).*(rhoP*hp/Qp/I)))/sqrt(2) ;
        plot3(freq,real(kbs2),imag(kbs2),'-','linewidth',1) ;
    % Torsion wave
        kt = omega.*(rhoP/Gp)^.5 ;
        plot3(freq,real(kt),imag(kt),'-','linewidth',1) ;
% Computed wavenumbers
    FREQ = freq + 0*K0 ;
    Kdir = K0(:).*dir ;
    plothandle = plot3(FREQ(:),real(K0(:)),imag(K0(:)),'.r','markersize',10,'linewidth',.1) ;
%     plot3(freq,real(K0),imag(K0),'+k','markersize',4,'linewidth',.1) ;
%     plot3(freq,real(Ks),imag(Ks),'xr','markersize',4,'linewidth',.1) ;
% Set plot options
    pkg.fem.bloch.setPlot(gammaMax,plotType,logScale) ;
%     pkg.fem.bloch.waveModeAnimation(mesh,Kdir,U0,plothandle) ;