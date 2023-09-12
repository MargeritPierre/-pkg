%% WAVES IN AN HOMOGENEOUS PRISMATIC WAVEGUIDE
% ESTIMATION OF A HYSTERETIC MODEL FROM SYNTHESIZED DATA
clc
close all
clear all
%% SYNTHESIZE THE DISPERSION CHARACTERISTICS
% Parameters
    % Geometry
        L = [4 .5] ; % unit cell size (mm)
        dx = .24*[4 1] ; % element size (mm)
    % Material
        E0 = 3.5e3*(1+.05i) ; % Young modulus (MPa)
        G0 = E0/2.6 ; % shear modulus (MPa)
        Kappa0 = E0.*G0./(3*(3*G0-E0)) ;
        rho = 1220e-12 ; % material density (tons/mm^3)
    % Wave 
        freq0 = logspace(4,6,100) ; % wave frequency
        nModes = 100 ; % number of extracted wave modes
        dir = [NaN NaN 1] ;
    % Display
        plotType = 'k' ; % 'c', 'k' , 'l' or 'g'
        logScale = true ;
        gammaMax = .1 ; inf ;
% Build the mesh
    mesh = pkg.geometry.mesh.GridMesh(L,dx) ;
    if numel(L)==2 ; mesh.setElementTypes(pkg.geometry.mesh.elements.Quad8) ; end
% Display the mesh
    clf ; axis equal tight ; view([30 30]) ;
    plot(mesh) ;
% Build the FEM matrices
    C0 = pkg.fem.bloch.stiffness(E0,G0) ;
    [K00,K0i,Kij,M,P] = pkg.fem.bloch.FEM(mesh,C0,rho,[],[]) ;
% Compute the wavenumbers and modes
    [k0,U0,omega] = pkg.fem.bloch.solve(K00,K0i,Kij,M,freq0,dir,nModes) ;
    [k0,U0] = pkg.fem.bloch.sort(k0,U0) ;
    U0 = reshape(P*U0(:,:),[mesh.nNodes 3 size(U0,2:ndims(U0))]) ;
%% DISPLAY RESULTS
    clf ;
% Computed wavenumbers
    FREQ = freq0 + 0*k0 ;
    Kdir = k0(:).*dir ;
    %plothandle = plot3(FREQ(:),real(k0(:)),imag(k0(:)),'ok','markersize',4,'linewidth',.1) ;
    plothandle = scatter3(FREQ(:),real(k0(:)),imag(k0(:)),25,'r','filled') ;
% Set plot options
    pkg.fem.bloch.setPlot(gammaMax,plotType,logScale) ;
    %pkg.fem.bloch.waveModeAnimation(mesh,Kdir,U0,plothandle,true)
    
%% ELEMENTARY MATRICES
Kappa = Kappa0*1.00001 ; G = G0*0.99999 ; 
% C = Kappa*Ck + G*Cg (Kappa & G are isostatic and shear moduli)
    Ck = blkdiag(ones(3),zeros(3)) ; 
    Cg = blkdiag(2*eye(3),eye(3)) - 2/3*Ck ;
% Elementary matrices in the third direction (wavenumber direction)
    [K00_k,K0i_k,Kij_k] = pkg.fem.bloch.FEM(mesh,Ck,rho,[],[]) ;
    K1_k = K0i_k{3} ; K2_k = Kij_k{3,3} ;
    [K00_g,K0i_g,Kij_g] = pkg.fem.bloch.FEM(mesh,Cg,rho,[],[]) ;
    K1_g = K0i_g{3} ; K2_g = Kij_g{3,3} ;
% Polynomial eigenvalue matrices
    K00 = Kappa*K00_k + G*K00_g ;
    K1 = Kappa*K1_k + G*K1_g ;
    K2 = Kappa*K2_k + G*K2_g ;
% Compute wavenumbers
    omega = 2*pi*freq0 ; nFreq = numel(freq0) ;
    nDOFs = size(K00,1) ; I = speye(nDOFs) ; O = sparse(nDOFs,nDOFs) ;
    k = NaN(nModes,nFreq) ;
    U = NaN(nDOFs,nModes,nFreq) ;
    wtbr = waitbar(0,'Bloch Wave Computation..') ;
    for ww = 1:nFreq
        K0 = omega(ww)^2*M-K00 ;
        A = blkdiag(K0,I) ; 
        B = [K1 K2 ; I O] ;
        [uu,ku] = eigs(A,B,nModes,'sm') ;
        k(:,ww) = diag(ku) ;
        U(:,:,ww) = reshape(uu(1:end/2,:),nDOFs,nModes) ;
        wtbr = waitbar(ww/nFreq,wtbr) ;
    end
    delete(wtbr) ;
    [k,U] = pkg.fem.bloch.sort(k,U) ;
    %U = reshape(P*U(:,:),[mesh.nNodes 3 size(U,2:ndims(U))]) ;
%% Wavenumber gradient dki/dp = -(Ui'*Ri*Ui)/(Ui'*Si*Ui)
% with Ri = (ki^2*dK2/dp + ki*dK1_dp + dK0_dp)
%  and Si = (2ki*K2 + K1)
    u = reshape(U,[],nModes*nFreq) ;
    uk = u.*k(:).' ;
    uk2 = u.*(k(:).^2).' ;
    w2Mu = M*(u.*repelem(omega.^2,50)) ;
    uSu = sum((2*K2*uk + K1*u).*conj(u),1) ;
    uRu_k = sum((K2_k*uk2 + K1_k*uk + K00_k*u).*conj(u),1) ;
    uRu_g = sum((K2_g*uk2 + K1_g*uk + K00_g*u).*conj(u),1) ;
    dk_dKappa = -reshape(uRu_k./uSu,nModes,nFreq) ;
    dk_dG = -reshape(uRu_g./uSu,nModes,nFreq) ;
    dk_dp = [dk_dKappa(:) dk_dG(:)] ;
    p = [Kappa ; G] ; p0 = [Kappa0 ; G0] ;
%% Cul invalid points
    valid = true(size(k(:))) ;
    valid = valid & (abs(k(:)-k0(:))./abs(k0(:)))<.0001 ;
    valid = valid & abs(imag(k(:)))./abs(real(k(:)))<1/10 ;
    
    dp = -(dk_dp(valid,:)'*dk_dp(valid,:))\(dk_dp(valid,:)'*(k(valid)-k0(valid))) ;
    [p0 p p0-p dp] 
%% Display
    displayTag = 'candidate' ; delete(findall(gcf,'tag',displayTag)) ; 
    FREQ = freq0 + 0*k ;
    Kdir = k(:).*dir ;
    plot3(FREQ(:),real(k(:)),imag(k(:)),'.r','markersize',8,'linewidth',.1,'tag',displayTag) ;
    %quiver3(FREQ(:),real(k(:)),imag(k(:)),FREQ(:)*0,real(dk_dKappa(:)),imag(dk_dKappa(:)),'autoscale',false,'tag',displayTag) ;
    pkg.fem.bloch.setPlot(gammaMax,plotType,logScale) ;


%% Thick plate reduced model
nu0 = E0./(2*G0)-1 ;
Q0 = E0/(1-nu0^2) ;
k = reshape(logspace(-1,1,1000),[],1) ;
w = 2*pi*reshape(logspace(log10(min(freq0)),log10(max(freq0)),1000),1,[]) ;
kl = sqrt(rho*w.^2/Q0 - k.^2) ;
kt = sqrt(rho*w.^2/G0 - k.^2) ;

ll = L(1)/2.*kl ; tt = L(1)/2.*kt ;
aa = 4*(k.^2.*kl.*kt) ; bb = (kt.^2-k.^2).^2 ;

detAs = aa.*sin(ll).*cos(tt) + bb.*sin(tt).*cos(ll) ;
detAs = abs(detAs)./medfilt2(abs(detAs),10*[1 1],'symmetric') ;

detAa = aa.*sin(tt).*cos(ll) + bb.*sin(ll).*cos(tt) ;
detAa = abs(detAa)./medfilt2(abs(detAa),10*[1 1],'symmetric') ;

xi = pi^2/12 ; eta = 12/(L(2)^2) ; 
ssss = sqrt((1/xi/G0 - 1/Q0)^2 + 4*eta/rho/Q0./(w.^2)) ;
mub2 = rho*w.^2/2.*((1/xi/G0 + 1/Q0) + ssss) ;
musa2 = rho*w.^2/2.*((1/xi/G0 + 1/Q0) - ssss) ;
must2 = rho*w.^2/G0-eta*xi ;
kb = sqrt(mub2-k.^2) ; cb = cos(L(1)/2.*kb) ; sb = sin(L(1)/2.*kb) ; 
ksa = sqrt(musa2-k.^2) ; csa = cos(L(1)/2.*ksa) ; ssa = sin(L(1)/2.*ksa) ;
kst = sqrt(must2-k.^2) ; cst = cos(L(1)/2.*kst) ; sst = sin(L(1)/2.*kst) ;

Sb3 = eta*xi*G0*((nu0-1)*k.^2 + mub2) ;
Ssa3 = eta*xi*G0*((nu0-1)*k.^2 + musa2) ;
Ab1 = kb.*(rho*w.^2-Q0*mub2) ;
Ab2 = 2*eta*xi*G0.*k.*kb ;
Asa1 = ksa.*(rho*w.^2-Q0*musa2) ;
Asa2 = 2*eta*xi*G0.*k.*ksa ;
Sst1 = k ;
Sst2 = k.^2 - kst.^2 ;
Ast3 = (1-nu0)*k.*kst ;


Ast3Ab1Asa2mAst3Ab2Asa1 = ...Ast3.*(Ab1.*Asa2-Ab2.*Asa1) ...
                          -4*G0*rho*w.^2.*k.^2.*kst.*kb.*ksa.*ssss ...
                          ; 
Ssa3Ab2Sst1mSsa3Ab1Sst2 = ...Ssa3.*(Ab2.*Sst1-Ab1.*Sst2) ...
                          kb.*(ksa.^2 + nu0*k.^2).*((rho*w.^2-2*G0*k.^2).*(k.^2+kst.^2) + Q0*(k.^2 + kb.^2).*(k.^2-kst.^2)) ...
                          ; 
Sb3Asa1Sst2mSb3Asa2Sst1 = ...Sb3.*(Asa1.*Sst2-Asa2.*Sst1) ...
                          -ksa.*(kb.^2 + nu0*k.^2).*((rho*w.^2-2*G0*k.^2).*(k.^2+kst.^2) + Q0*(k.^2 + ksa.^2).*(k.^2-kst.^2)) ...
                          ...-ksa.*(mub2-2*G0/Q0*k.^2).*(2*eta*xi*G0*k.^2 + (Q0*musa2-rho*w.^2).*(k.^2-kst.^2)) ...
                          ; 
                      
detBs = cst.*sb.*ssa.*Ast3Ab1Asa2mAst3Ab2Asa1 ...
      + csa.*sb.*sst.*Ssa3Ab2Sst1mSsa3Ab1Sst2 ...
      + cb.*ssa.*sst.*Sb3Asa1Sst2mSb3Asa2Sst1 ;
detBs = abs(detBs)./medfilt2(abs(detBs),floor(size(detBs)/100),'symmetric') ;

detBa = sst.*cb.*csa.*Ast3Ab1Asa2mAst3Ab2Asa1 ...
      + ssa.*cb.*cst.*Ssa3Ab2Sst1mSsa3Ab1Sst2 ...
      + sb.*csa.*cst.*Sb3Asa1Sst2mSb3Asa2Sst1 ;
detBa = abs(detBa)./medfilt2(abs(detBa),floor(size(detBa)/100),'symmetric') ;


clf reset ; 
plot3(FREQ(:),real(k0(:)),imag(k0(:)),'.r','markersize',3,'linewidth',.1) ;
% plot3(w/2/pi,real(sqrt(mub2)),eps+imag(kb),'ok','markersize',5,'linewidth',.1) ;
% plot3(w/2/pi,real(sqrt(musa2)),eps+imag(kb),'ok','markersize',5,'linewidth',.1) ;
% plot3(w/2/pi,real(sqrt(must2)),eps+imag(kb),'ok','markersize',5,'linewidth',.1) ;
pkg.fem.bloch.setPlot(gammaMax,plotType,logScale) ;
srf = surf(repmat(w/2/pi,[numel(k) 1])...
            ,repmat(k,[1 numel(w)])...
            ,eps+0*w.*k...
            ,(abs(  detBa.*detBs.*detAa.*detAs )).*ones(size(detAs))) ;
srf.FaceColor = 'interp' ; srf.EdgeColor = 'none' ;
% axis tight
set(gca,'colorscale','log')
caxis([.5 1]) ;
colormap((gray))
set(gca,'xscale','log','yscale','log')














    