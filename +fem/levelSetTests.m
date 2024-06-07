%% SOME TESTS ON LEVELSET METHODS

%% DISTANCE FUNCTION, FINITE DIFFERENCES, NORMAL EVOLUTION
% see Dapogny & Maitre, "An introduciton to the level-set method"
% dphi_dt + V.grad(phi) = 0
% V = v*N normal velocity
% N = grad(phi)/|grad(phi)| outgoing normal
import pkg.geometry.levelset.*

% Parameters
D = 1 ;
R = D/4 ;
H = R/2 ;
m = R ;
dx = R/30 ;
dt = .5*dx^2 ;
T = .02 ;
nReDistIt = 1 ;
dtDist = 50*dt ; 0.001 ;
epsDist = 3*dx ;
plotFreq = 10 ;
nIt = round(T/dt) ;

lvlst = Circle([0 0],R) + Circle([D 0],R) + Rectangle([0 -H/2 ; D H/2]) ;
%lvlst = Rectangle([0 0 ; D H]) ;

% Derivatives: centered, after, before
dCen = @(f)reshape([f(2,:)-f(1,:) ; .5*(f(3:end,:)-f(1:end-2,:)) ; f(end,:)-f(end-1,:)],size(f)) ;
dPlus = @(f)reshape(f([2:end,end],:)-f([1:end-1,end-1],:),size(f)) ;
dMinus = @(f)reshape(f([2,2:end],:)-f([1,1:end-1],:),size(f)) ;
diff = @(f,dim,der) permute(der(permute(f,[dim 1:dim-1 dim+1:ndims(f)])),[2:dim 1 dim+1:ndims(f)])/dx ;


bbox = lvlst.BoundingBox+[-1;1]*m ;
[XX,YY] = meshgrid(bbox(1,1):dx:bbox(2,1),bbox(1,2):dx:bbox(2,2)) ;
PHI = reshape(lvlst.Function([XX(:) YY(:)]),size(XX)) ;

cla ; axis equal
srf = surf(XX,YY,PHI,'Facecolor','interp','EdgeColor','none') ;
[~,ctr] = contour(XX,YY,PHI,'LevelList',0,'linecolor','k','linewidth',2) ;
colorbar
drawnow

plotTime = tic ;
for it = 1:nIt
% Velocity field
    dPHI_dx = diff(PHI,2,dCen) ;
    dPHI_dy = diff(PHI,1,dCen) ;
    gradPHI = cat(3,dPHI_dx,dPHI_dy) ;
    normGradPHI = sqrt(sum(gradPHI.^2,3)) ;
    normal = gradPHI./(normGradPHI+eps) ;
    curv = diff(normal(:,:,1),2,dCen) + diff(normal(:,:,2),1,dCen) ;
    V = -(curv) ;
% Level-set advection
    Dpx = diff(PHI,2,dPlus) ;
    Dmx = diff(PHI,2,dMinus) ;
    Dpy = diff(PHI,1,dPlus) ;
    Dmy = diff(PHI,1,dMinus) ;
    Gp = sqrt(max(max(Dmx,0),-min(Dpx,0)).^2+max(max(Dmy,0),-min(Dpy,0)).^2) ;
    Gm = sqrt(max(max(Dpx,0),-min(Dmx,0)).^2+max(max(Dpy,0),-min(Dmy,0)).^2) ;
    PHI = PHI - sum(max(V,0).*Gp + min(V,0).*Gm,3)*dt ;
% Level-set reinitialization        
    signPHI = PHI./sqrt(PHI.^2+epsDist.^2) ;
    for dit = 1:nReDistIt
        Dpx = diff(PHI,2,dPlus) ;
        Dmx = diff(PHI,2,dMinus) ;
        Dpy = diff(PHI,1,dPlus) ;
        Dmy = diff(PHI,1,dMinus) ;
        Gp = sqrt(max(max(Dmx,0),-min(Dpx,0)).^2+max(max(Dmy,0),-min(Dpy,0)).^2) ;
        Gm = sqrt(max(max(Dpx,0),-min(Dmx,0)).^2+max(max(Dpy,0),-min(Dmy,0)).^2) ;
        PHI = PHI - dtDist*sum(max(signPHI,0).*(Gp-1) + min(signPHI,0).*(Gm-1),3) ;
    end
% Display
    if it==nIt || toc(plotTime)>1/plotFreq
        srf.ZData = PHI ;
        ctr.ZData = PHI ;
        drawnow ;
        plotTime=tic ;
    end
end



%% CHARACTERISTIC FUNCTION, FINITE DIFFERENCES
% see Olsson & Kreiss, "A conservative level set method for two phase flow"
% phi = 1/(1+exp(d/eps)) thick heaviside function (levelset at 0.5)
% dphi_dt + V.grad(phi) = 0 advection equation with V velocity vector
% dphi_dtau + div(phi.*(1-phi).*N) = eps*div(grad(phi)) reinitialization
% with N = grad(phi)/|grad(phi)| the ingoing normal
import pkg.geometry.levelset.*

% Parameters
D = 1 ;
R = D/4 ;
H = R/1 ;
m = 2.5*R ;
dx = R/20 ;
dt = 5*dx^2 ;
T = 1 ; %.02 ;
nReDistIt = 5 ;
dd = 0.1 ; beta = 1/2 ;
dtDist = beta*dx^(1+dd) ; 50*dt ; 
epsDist = beta*dx^(1-dd) ; 1*dx ; sqrt(2)*dx ;
plotFreq = 15 ;
nIt = round(T/dt) ;
projection = 'normal' ; % 'full' or 'normal'

lvlst = Circle([0 0],R) + Circle([D 0],R) + Rectangle([0 -H/2 ; D H/2]) ;
%lvlst = Rectangle([0 0 ; D H]) ;

% Derivatives: centered, after, before
dCen = @(f)reshape([f(2,:)-f(1,:) ; .5*(f(3:end,:)-f(1:end-2,:)) ; f(end,:)-f(end-1,:)],size(f)) ;
dPlus = @(f)reshape(f([2:end,end],:)-f([1:end-1,end-1],:),size(f)) ;
dMinus = @(f)reshape(f([2,2:end],:)-f([1,1:end-1],:),size(f)) ;
diff = @(f,dim,der) permute(der(permute(f,[dim 1:dim-1 dim+1:ndims(f)])),[2:dim 1 dim+1:ndims(f)])/dx ;


bbox = lvlst.BoundingBox+[-1;1]*m ;
[XX,YY] = meshgrid(bbox(1,1):dx:bbox(2,1),bbox(1,2):dx:bbox(2,2)) ;
PHI = reshape(lvlst.Function([XX(:) YY(:)]),size(XX)) ;

PHI = 1./(1+exp(PHI./epsDist)) ;

cla ; axis equal
srf = surf(XX,YY,PHI,'Facecolor','interp','EdgeColor','none') ;
[~,ctr] = contour3(XX,YY,PHI,'LevelList',0.5,'linecolor','k','linewidth',2) ;
colorbar
drawnow ;

plotTime = tic ;
for it = 1:nIt
% Velocity field
    dPHI_dx = diff(PHI,2,dCen) ;
    dPHI_dy = diff(PHI,1,dCen) ;
    gradPHI = cat(3,dPHI_dx,dPHI_dy) ;
    normGradPHI = sqrt(sum(gradPHI.^2,3)) ;
    normal = gradPHI./(normGradPHI+eps) ;
    curv = diff(normal(:,:,1),2,dCen) + diff(normal(:,:,2),1,dCen) ;
    V = -curv.*normal ;
    V = 2*pi*cat(3,-(YY-mean(YY(:))),(XX-mean(XX(:)))) ;
% Level-set advection
    switch projection
        case 'normal' 
            V = sum(V.*normal,3) ;
            Dpx = diff(PHI,2,dPlus) ;
            Dmx = diff(PHI,2,dMinus) ;
            Dpy = diff(PHI,1,dPlus) ;
            Dmy = diff(PHI,1,dMinus) ;
            Gp = sqrt(max(max(Dmx,0),-min(Dpx,0)).^2+max(max(Dmy,0),-min(Dpy,0)).^2) ;
            Gm = sqrt(max(max(Dpx,0),-min(Dmx,0)).^2+max(max(Dpy,0),-min(Dmy,0)).^2) ;
            dPHI_dt = - sum(max(V,0).*Gp + min(V,0).*Gm,3) ;
        case 'full' 
            Fp = .5*(PHI(:,[2:end,end])+PHI).*V(:,:,1) ;
            Fm = .5*(PHI(:,[1,1:end-1])+PHI).*V(:,:,1) ;
            Gp = .5*(PHI([2:end,end],:)+PHI).*V(:,:,2) ;
            Gm = .5*(PHI([1,1:end-1],:)+PHI).*V(:,:,2) ;
            dPHI_dt = (-1/dx)*(Fp-Fm+Gp-Gm) ;
    end
    dPHI_dt([1,end],:) = 0 ; dPHI_dt(:,[1,end]) = 0 ; 
    PHI = PHI + dPHI_dt*dt ;
% Level-set reinitialization   
    dPHI_dx = diff(PHI,2,dCen) ;
    dPHI_dy = diff(PHI,1,dCen) ;
    gradPHI = cat(3,dPHI_dx,dPHI_dy) ;
    normGradPHI = sqrt(sum(gradPHI.^2,3)) ;
    normal = gradPHI./(normGradPHI+eps) ;     
    for dit = 1:nReDistIt
        f = PHI.*(1-PHI).*normal ;
        Fp = .5*(f(:,[2:end,end],1)+f(:,:,1))-(epsDist/dx)*(PHI(:,[2:end,end])-PHI(:,[1:end-1,end-1])) ;
        Fm = .5*(f(:,[1,1:end-1],1)+f(:,:,1))-(epsDist/dx)*(PHI(:,[2,2:end])-PHI(:,[1,1:end-1])) ;
        Gp = .5*(f([2:end,end],:,2)+f(:,:,2))-(epsDist/dx)*(PHI([2:end,end],:)-PHI([1:end-1,end-1],:)) ;
        Gm = .5*(f([1,1:end-1],:,2)+f(:,:,2))-(epsDist/dx)*(PHI([2,2:end],:)-PHI([1,1:end-1],:)) ;
        dPHI_dTau = (-1/dx)*(Fp-Fm+Gp-Gm) ;
        dPHI_dTau([1,end],:) = 0 ; dPHI_dTau(:,[1,end]) = 0 ; 
        PHI = PHI + dtDist*dPHI_dTau ;
    end
% Display
    if it==nIt || toc(plotTime)>1/plotFreq
        srf.ZData = PHI ;
        ctr.ZData = PHI ;
        drawnow ;
        plotTime = tic ;
    end
end



%% CHARACTERISTIC FUNCTION, VARIATIONAL FORMULATION
% see Olsson & Kreiss, "A conservative level set method for two phase flow"
% phi = 1/(1+exp(d/eps)) thick heaviside function (levelset at 0.5)
% 1) dphi_dt + V.grad(phi) = 0 advection equation with V velocity vector
% 2) dphi_dtau + div(phi.*(1-phi).*N) = eps*div(grad(phi)) reinitialization
% with N = grad(phi)/|grad(phi)| the INgoing normal
% Variationnal formulation: with "s" virtual field, assume phi=0 on boundaries
% 1) advection: int(s*dphi_dt*dx) - int((grad(s).u)*phi*dx) = 0 
%       assumes u.n = 0 on the domain boundaries (free slip BC) 
%       or vanishing phi on BCs
import pkg.geometry.levelset.*

% Parameters
D = 1 ;
R = D/4 ;
H = R/1 ;
m = 2.5*R ;
dx = R/10 ;
dt = .25*dx ; .5*dx^2 ;
T = 1 ; .02 ;
nReDistIt = 1 ;
dd = 0.1 ; beta = 1/2 ;
dtDist = beta*dx^(1+dd) ; 50*dt ; 
epsDist = beta*dx^(1-dd) ; 1*dx ; sqrt(2)*dx ;
plotFreq = 1 ;
nIt = round(T/dt) ;
theta = .5 ; % theta-scheme parameter: theta==1 -> implicit, theta=0 -> explicit

lvlst = Circle([0 0],R) + Circle([D 0],R) + Rectangle([0 -H/2 ; D H/2]) ;
%lvlst = Rectangle([0 0 ; D H]) ;

bbox = lvlst.BoundingBox+[-1;1]*m ;
[XX,YY] = meshgrid(bbox(1,1):dx:bbox(2,1),bbox(1,2):dx:bbox(2,2)) ;

mesh = pkg.geometry.mesh.Mesh(cat(3,XX,YY)) ;

PHI = reshape(lvlst.Function([XX(:) YY(:)]),size(XX)) ;
PHI = 1./(1+exp(PHI./epsDist)) ;

[ee,w,ie] = mesh.integration() ;
nGP = numel(w) ; W = diag(sparse(w)) ;
N = mesh.interpMat(ee,ie) ;
G = mesh.diffMat(ee,ie) ;
I = N'*W*N ;
I = diag(sum(I,2)) ; % mass lumping..
GG = cat(1,G{:})'*kron(speye(2),W)*cat(1,G{:}) ; % grad(s)*grad(phi)

keepDOF = ~mesh.boundaryNodes ; % phi=0 on boundaries..

cla ; axis equal
srf = plot(mesh,'EdgeColor','none','Deformation',[0 0 1].*PHI(:)) ;
[~,ctr] = contour3(XX,YY,PHI,'LevelList',0.5,'linecolor','k','linewidth',2) ;
colorbar
drawnow ;

plotTime = tic ;
for it = 1:nIt
% Velocity field
    V = pi*[-(YY(:)-mean(YY(:))) (XX(:)-mean(XX(:)))] ;
% Level-set advection
    u = N*sparse(V) ;
    gSu = diag(u(:,1))*G{1} + diag(u(:,2))*G{2} ;
    A = (1/dt)*I - theta*gSu'*W*N ;
    b = (I*PHI(:))/dt + (1-theta)*(gSu'*(W*(N*PHI(:)))) ;
    PHI(keepDOF) = A(keepDOF,keepDOF)\b(keepDOF) ;
% Level-set reinitialization   
    gradPHI = reshape(cat(1,G{:})*PHI(:),nGP,2) ;
    normGradPHI = sqrt(sum(gradPHI.^2,2)) ;
    normal = gradPHI./(normGradPHI+eps) ;  
    gSn = W*[diag(sparse(normal(:,1))) diag(sparse(normal(:,2)))]*cat(1,G{:}) ;
    for dit = 1:nReDistIt
        A = (1/dt)*I - theta*gSn'*(N.*(1-2*PHI(:))') + epsDist*GG ;
        b = gSn'*(N*(PHI(:).*(1-theta-(1-2*theta)*PHI(:)))) - epsDist*(1-theta)*GG*PHI(:) + (I*PHI(:))/dt ;
        PHI(keepDOF) = A(keepDOF,keepDOF)\b(keepDOF) ;
    end 
% Display
    if it==nIt || toc(plotTime)>1/plotFreq
        srf.Deformation = [0 0 1].*PHI(:) ;
        srf.CData = PHI(:) ;
        ctr.ZData = PHI ;
        drawnow ;
        plotTime = tic ;
    end
end
