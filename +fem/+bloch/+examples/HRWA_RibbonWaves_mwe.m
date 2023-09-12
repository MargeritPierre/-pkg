%% ESTIMATE MATERIAL PARAMETERS FROM HRWA RESULTS
% Assumes the waveguide to be a ribbon (bounded thick plate)

%% DATA LOADING USING A HRWA RESULT FIGURE
% to delete data, using brushing then "remove brushed" (press delete key)
figK = open('miniHP_K.fig') ;
axK = figK.Children(1) ;

%% EXTRACT THE VALID DATA
sca = findall(axK,'type','scatter') ;
sca.SizeData = 5 ; 
DATA = get(sca,{'xdata','ydata','zdata'}) ;
DATA = cat(1,DATA{:}).' ; % [F ReK ImK]

valid = all(~isnan(DATA),2) ;

F = DATA(valid,1) ;
K = (DATA(valid,2)+1i*DATA(valid,3))/1000 ;

%% INITIALIZE THE MODEL
% Ribbon geometry
b = 372/100 ; % width (mm)
h = 48/100 ; % height (mm)
% Material
E = 38e2*(1+0.00i) ; % Young modulus (MPa)
G = 14e2*(1+0.00i) ; E/2/(1+36/100) ; % shear modulus (MPa)
rho = 1220e-12 ; % material density (tons/mm^3)
nu = E./(2*G)-1 ;
Q = E/(1-nu^2) ;
% Determinant options
normalize = true ;
delta = 1e-4 ; % finite difference relative step for derivative estimation
dQ = delta*real(Q) ;
dG = delta*real(G) ;
dK = delta*abs(K) ;
% Display options
Npix = 500 ;
filtSz = 5.*[1 1] ;
wplot = 2*pi*logspace(log10(min(F)),log10(max(F)),Npix) ; % frequency
kplot = reshape(logspace(log10(min(real(K))),log10(max(real(K))),Npix),[],1) ; % wavenumber

Dplot = determinant(h,b,Q,G,rho,wplot,kplot,normalize) ;
Dplot = Dplot./medfilt2(abs(Dplot),filtSz,'symmetric') ;

tag = 'detplot' ;
delete(findall(figK,'tag',tag)) ;
srf = surf(axK ...
            ,repmat(wplot/2/pi,[numel(kplot) 1])...
            ,1000*repmat(kplot,[1 numel(wplot)])...
            ,min(get(gca,'zlim'))+0*log10(abs(Dplot)) ...
            ,(abs(Dplot)).*ones(size(Dplot))...
            ,'edgecolor','none' ...
            ,'facecolor','interp' ...
            ,'tag',tag ...
        ) ;
set(gca,'colorscale','log')
caxis([min(abs(srf.CData(:))) max(abs(srf.CData(:)))])
caxis([7e-1 1])
colormap((1-flip(gray,1).*[1 0 1]))

%%
for it = 1:1
    
%% Determinant values & derivatives at experimantal points
Dexp = determinant(h,b,Q,G,rho,2*pi*F,K,normalize) ;
dD_dQ = (determinant(h,b,Q+dQ,G,rho,2*pi*F,K,normalize) - determinant(h,b,Q-dQ,G,rho,2*pi*F,K,normalize))/2/dQ ;
dD_dG = (determinant(h,b,Q,G+dG,rho,2*pi*F,K,normalize) - determinant(h,b,Q,G-dG,rho,2*pi*F,K,normalize))/2/dG ;
 
    
% Optimize the distance (K0-K).^2 ?
if 1
    % Find the closest zeros of the determinant
    K0 = K ;
    for ii = 1:10 % descent algorithm
        D0 = determinant(h,b,Q,G,rho,2*pi*F,K0,normalize) ;
    % linear continuation D(K+dK) ~ D(K) + dK*D' = 0
        dD_dK0 = (determinant(h,b,Q,G,rho,2*pi*F,K0+dK,normalize) ...
                    - determinant(h,b,Q,G,rho,2*pi*F,K0-dK,normalize))./2./dK ;
        dK1 = -D0./dD_dK0 ;
    % Quadratic continuation D(K+dK) ~ D(K) + dK*D' + .5*dK^2*D'' = 0
        d2D_dK02 = (determinant(h,b,Q,G,rho,2*pi*F,K0+dK,normalize) ...
                    - 2*D0 + determinant(h,b,Q,G,rho,2*pi*F,K0-dK,normalize))./(dK.^2) ;
        aa = 2*D0./d2D_dK02 ; bb = dD_dK0./d2D_dK02 ; dd = sqrt(bb.^2-aa) ;
        dK2a = -bb + dd ;
        dK2b = -bb - dd ;
    % Take the smallest step
        dK0 = dK2a.*(abs(dK2b)>abs(dK2a)) + dK2b.*(abs(dK2b)<=abs(dK2a)) ;
        K0 = K0 + 1*dK0 ;
    end
    tag = 'nearestK' ;
    delete(findall(figK,'tag',tag)) ;
    plot3(findobj(figK,'type','axes')...
            ,[F F].',real([K K0].')*1000,imag([K K0].')*1000 ...
            ,'o-m' ...
            ,'linewidth',.1 ...
            ,'markersize',3 ...
            ,'tag',tag) ;
   
    % Determinant gradient at the solution
    dD_dQ = (determinant(h,b,Q+dQ,G,rho,2*pi*F,K0,normalize) - determinant(h,b,Q-dQ,G,rho,2*pi*F,K0,normalize))/2/dQ ;
    dD_dG = (determinant(h,b,Q,G+dG,rho,2*pi*F,K0,normalize) - determinant(h,b,Q,G-dG,rho,2*pi*F,K0,normalize))/2/dG ;
    dD_dK0 = (determinant(h,b,Q,G,rho,2*pi*F,K0+dK,normalize) ...
                - determinant(h,b,Q,G,rho,2*pi*F,K0-dK,normalize))./2./dK ;
    dD_dQ = -dD_dQ./dD_dK0 ;
    dD_dG = -dD_dG./dD_dK0 ;
    Dexp = (K0-K) ;
end

% weighting
%weights = abs(K).^2.*F.^-4.*abs(Dexp).^-1.5 ;
weights = abs(K).^0.*F.^-2.*abs(Dexp).^-0 ;
W = diag(sparse(weights)) ;

% frequency dependence
switch 1
    case 0 % no frequency dependence, hysteretic Q and G
        Q = mean(Q(:)) ; G = mean(G(:)) ;
        dD_dP = [dD_dQ(:) dD_dG(:)] ;
        T = 1 ;
    case 1 % regularized frequency dependence
        % all parameters independent
        Q = Q + 0*F ; G = G+0*F ;
        dD_dP = [diag(sparse(dD_dQ(:))) diag(sparse(dD_dG(:)))] ;
        % one parameter by unique frequency
        [uF,ia,ic] = unique(F(:),'sorted') ;
        T = sparse(1:numel(F),ic(:)',1,numel(F),numel(uF)) ;
        dD_dP = dD_dP*blkdiag(T,T) ;
        % frequency regularization
        beta = 1e3 ;
        % ddw2*f = d2f_dw2
        ss = speye(numel(uF)-2,numel(uF)) ;
        ddw2 = 2*( (circshift(ss,2,2)-circshift(ss,1,2))./(uF(3:end)-uF(2:end-1)) - (circshift(ss,1,2)-ss)./(uF(2:end-1)-uF(1:end-2)) ) ./ (uF(3:end)-uF(1:end-2)); 
        % associated jacobian & cost
        jreg = blkdiag(ddw2,ddw2) ;
        dD_dP = [dD_dP ; jreg] ;
        Dexp = [Dexp(:) ; jreg*[Q(ia);G(ia)]] ;
        % weights
        Wd2w = beta*speye(2*(numel(uF)-2)) ;
        W = blkdiag(W,Wd2w) ;
end

% solve
dP = -.5*((dD_dP'*W*dD_dP)\(dD_dP'*W*Dexp(:))) ;
P = [Q;G] + blkdiag(T,T)*dP ;

% contraints
if isreal(K)
    P = abs(real(P)).*(1+0.000000001i) ;
else
    P = abs(real(P))+1i*abs(imag(P)) ;
end

% apply
Q = P(1:end/2) ; G = P(end/2+1:end) ;

% secondary results
nu = 1-2*G./Q ;
E = Q.*(1-nu.^2) ;


Dplot = determinant(h,b,mean(Q(:)),mean(G(:)),rho,wplot,kplot,normalize) ;
Dplot = Dplot./medfilt2(abs(Dplot),filtSz,'symmetric') ;

tag = 'detplot' ;
delete(findall(figK,'tag',tag)) ;
srf = surf(axK,repmat(wplot/2/pi,[numel(kplot) 1])...
            ,1000*repmat(kplot,[1 numel(wplot)])...
            ,min(get(gca,'zlim'))+0*log10(abs(Dplot)) ...
            ,(abs(Dplot)).*ones(size(Dplot))...
            ,'edgecolor','none' ...
            ,'facecolor','interp' ...
            ,'tag',tag ...
        ) ;
set(gca,'colorscale','log')
caxis([min(abs(srf.CData(:))) max(abs(srf.CData(:)))])
caxis([7e-1 1])
colormap((1-flip(gray,1).*[1 0 1]))
colorbar

tag = 'matplot' ;
figMat = findall(0,'tag',tag) ;
if isempty(figMat) 
    figMat = figure('tag',tag) ;
    axMat = gobjects(0) ;
    axMat(1) = axes('outerposition',[0 .5 1 .5]) ; xlabel Frequency ; ylabel 'Storage Moduli' ;
    axMat(2) = axes('outerposition',[0 0 1 .5]) ; xlabel Frequency ; ylabel 'Loss Factor' ; 
end
cla(axMat(1)) ;
    plot(axMat(1),F,real(E).*ones(size(F)),'.r') ;
    plot(axMat(1),F,real(G).*ones(size(F)),'.b') ;
cla(axMat(2)) ;
    plot(axMat(2),F,imag(E)./real(E).*ones(size(F)),'.r') ;
    plot(axMat(2),F,imag(G)./real(G).*ones(size(F)),'.b') ;

drawnow ;
end

%% DETERMINANT COMPUTATION
function [D,detAs,detAa,detBs,detBa] = determinant(h,b,Q,G,rho,w,k,normalize)
if nargin<8 || isempty(normalize) ; normalize = false ; end
nu = 1-2.*G./Q ;
% IN-PLANE MOTION
% wavenumbers
kl = sqrt(rho.*w.^2./Q - k.^2) ; % longitudinal
kt = sqrt(rho.*w.^2./G - k.^2) ; % transverse
% precompute
cl = cos(b./2.*kl) ; sl = sin(b./2.*kl) ; 
ct = cos(b./2.*kt) ; st = sin(b./2.*kt) ; 
Al = 4*(k.^2.*kl.*kt) ; At = (kt.^2-k.^2).^2 ;
% determinants
detAs = Al.*sl.*ct + At.*st.*cl ; % symmetrical modes
detAa = Al.*st.*cl + At.*sl.*ct ; % antisymmetrical modes
% normalized
if normalize
    detAs = detAs./Al...
            .*abs(k).^2.*abs(w).^-2 ...
            .*exp(b/2*imag(kt+kl)) ; % symmetrical modes
    detAa = detAa./Al...
            .*abs(k).^2.*abs(w).^-2 ...
            .*abs(kt).^1.*max(abs(k),abs(sqrt(rho./G).*w)).^-1 ...
            ....*kl.^-1.*max(k,sqrt(rho/Q).*w).^1 ...
            .*exp(b/2*imag(kt)).*exp(b/2*imag(kl)) ...
            ; % antisymmetrical modes
end

% OUT-OF-PLANE MOTION
xi = sqrt(pi^2/12) ; % shear correction factor
eta = 12./(h.^2) ; % shape factor
% wavenumbers (isotropic part mui2 = ki2 + k2)
ssss = sqrt((1./xi./G - 1./Q).^2 + 4*eta./rho./Q./(w.^2)) ; 
mub2 = rho.*w.^2./2.*((1./xi./G + 1./Q) + ssss) ; % bending
musa2 = rho.*w.^2./2.*((1./xi./G + 1./Q) - ssss) ; % axial shear
must2 = rho.*w.^2./G-eta.*xi ; % transverse shear
% width-directio wavenumbers & trigo
kb = sqrt(mub2-k.^2) ; cb = cos(b./2.*kb) ; sb = sin(b./2.*kb) ; 
ksa = sqrt(musa2-k.^2) ; csa = cos(b./2.*ksa) ; ssa = sin(b./2.*ksa) ;
kst = sqrt(must2-k.^2) ; cst = cos(b./2.*kst) ; sst = sin(b./2.*kst) ;
% components of the system matrix
% Sb3 = eta.*xi.*G.*((nu-1)*k.^2 + mub2) ;
% Ssa3 = eta.*xi.*G.*((nu-1)*k.^2 + musa2) ;
% Ab1 = kb.*(rho.*w.^2-Q.*mub2) ;
% Ab2 = 2.*eta.*xi.*G.*k.*kb ;
% Asa1 = ksa.*(rho.*w.^2-Q.*musa2) ;
% Asa2 = 2.*eta.*xi.*G.*k.*ksa ;
% Sst1 = k ;
% Sst2 = k.^2 - kst.^2 ;
% Ast3 = (1-nu).*k.*kst ;
% precomputation
Ast3Ab1Asa2mAst3Ab2Asa1 = ...Ast3.*(Ab1.*Asa2-Ab2.*Asa1) ...
                          -4.*G.*rho.*w.^2.*k.^2.*kst.*kb.*ksa.*ssss ...
                          ; 
Ssa3Ab2Sst1mSsa3Ab1Sst2 = ...Ssa3.*(Ab2.*Sst1-Ab1.*Sst2) ...
                          kb.*(ksa.^2 + nu.*k.^2).*((rho.*w.^2-2.*G.*k.^2).*(k.^2+kst.^2) + Q.*(k.^2 + kb.^2).*(k.^2-kst.^2)) ...
                          ; 
Sb3Asa1Sst2mSb3Asa2Sst1 = ...Sb3.*(Asa1.*Sst2-Asa2.*Sst1) ...
                          -ksa.*(kb.^2 + nu.*k.^2).*((rho.*w.^2-2.*G.*k.^2).*(k.^2+kst.^2) + Q.*(k.^2 + ksa.^2).*(k.^2-kst.^2)) ...
                          ; 
% determinants
detBs = cst.*sb.*ssa.*Ast3Ab1Asa2mAst3Ab2Asa1 ...
      + csa.*sb.*sst.*Ssa3Ab2Sst1mSsa3Ab1Sst2 ...
      + cb.*ssa.*sst.*Sb3Asa1Sst2mSb3Asa2Sst1 ; % symmetric modes

detBa = sst.*cb.*csa.*Ast3Ab1Asa2mAst3Ab2Asa1 ...
      + ssa.*cb.*cst.*Ssa3Ab2Sst1mSsa3Ab1Sst2 ...
      + sb.*csa.*cst.*Sb3Asa1Sst2mSb3Asa2Sst1 ; % antisymmetric modes
% normalized
if normalize
    detBs = detBs...
            ./abs(w.^3.*k.^0.*kst.^0.*kb.^1.*ksa.^0.*ssss.^0)...
            ./abs(k.^2+kst.^2) ...
            .*exp(b/2*(imag(kb+kst)-abs(imag(ksa)))) ...
            ; % symmetrical modes
    detBa = detBa ...
            ./abs(w.^3.*k.^0.*kst.^0.*kb.^0.*ksa.^0.*ssss.^0)...
            ./(k.^2+kst.^2) ...
            .*exp(b/2*(imag(kb+kst)-abs(imag(ksa)))) ...
            ; % antisymmetrical modes
end

% complete determinant
D = detAs.*detAa.*detBa.*detBs ;

end