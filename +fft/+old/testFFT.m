%% FFT SCHEME, CPU, N-gridded variables are of size [N , ...]
clearvars
profile on

L = [1 1] ;
N = 128*ones(size(L)) ;
D = numel(N) ;

dx = L./N ; 
dxi = 2*pi./L ;

idx = arrayfun(@colon,0*N,N-1,'UniformOutput',false) ;
[idx{:}] = ndgrid(idx{:}) ;
idx = cat(D+1,idx{:}) ;
% idx = reshape(idx,[],D) ;

x = idx.*reshape(dx,[ones(1,D) D]) ;
xi = (idx-reshape(floor(N/2),[ones(1,D) D])).*reshape(dxi,[ones(1,D) D]) ;
xi = shift_fft(xi,N,true) ;

x0 = .5*L ; R0 = .3*min(L) ;
xS = [0*L;L;eye(D).*L;(1-eye(D)).*L] ; Rs = R0 ;
contrast = 1e3 ;

r20 = sum((x-reshape(x0',[ones(1,D) D size(x0,1)])).^2,D+1) ;
phi0 = min(r20-reshape(R0.^2,[ones(1,D+1) numel(R0)]),[],D+2) ;
r2s = sum((x-reshape(xS',[ones(1,D) D size(xS,1)])).^2,D+1) ;
phiS = min(r2s-reshape(Rs.^2,[ones(1,D+1) numel(Rs)]),[],D+2) ;
phi = .5*(sign(phi0)-sign(phiS)) ;

kappa0 = 1 ; mu0 = 1 ; rho0 = 1 ;
modulation = sqrt(contrast).^phi ;
kappa = kappa0*modulation ;
mu = mu0*modulation ;
rho = rho0*modulation ;

C = assemble_isoC(kappa,mu) ;
C0 = median(C,[1:D]) ; %assemble_isoC(kappa0,mu0,D) ;

method = 'analytic' ;

Id = reshape(eye(D),[ones(1,D) D D]) ;
tic ; P0 = div(tensprod(C0,grad(Id,xi,method),2),xi,method) ; toc
singular = all(abs(P0)<sqrt(eps),D+(1:2)) ;
P0 = P0 + singular.*reshape(eye(D),[ones(1,D) D D]) ;
tic ; iP0 = permute(pageinv(permute(P0,[D+(1:D) 1:D])),[D+(1:D) 1:D]) ; toc

sig_op = @(eps)tensprod(C,eps,2) ;
sig_op = @(eps)apply_isoC(eps,kappa,mu) ;
A_op = @(u)div(to_fourier(sig_op(to_spatial(grad(u,xi,method),N)),N),xi,method) ;
iA0_op = @(u)tensprod(iP0,u,1) ;

A_op_dof = @(u)grid2dof(A_op(dof2grid(u,N,D))) ;
iA0_op_dof = @(u)grid2dof(iA0_op(dof2grid(u,N,D))) ;

switch D
    case 2
        E = .1*reshape([0 1 ; 1 0],[ones(1,D) D D]) ;
    case 3
        E = .1*reshape([0 1 0 ; 1 0 0 ; 0 0 0],[ones(1,D) D D]) ;
end
S = tensprod(C,E,2) ;
f = -div(to_fourier(S,N),xi,method) ;

tic ; u = pcg(A_op_dof,f(:),1e-3,100,iA0_op_dof) ; toc
u = dof2grid(u,N,D) ;
u = to_spatial(u,N) ;
eps = sym(to_spatial(grad(to_fourier(u,N),xi,method),N)) + E ;
u = u + tensprod(E,x,1) ;

profile off

switch D
    case 2
        clf ; axis equal tight
        clr = kappa ;
        imagesc(x(:,1,1),x(1,:,2),clr')
        %surf(x(:,:,1),x(:,:,2),clr')
        set(gca,'colorscale','log')
        colorbar
        clf ; axis equal tight off
        xp = x + real(u) ;
        clr = sqrt(sum(abs(eps).^2,[3 4])) ; sqrt(sum(abs(u).^2,3)) ; sqrt(sum(S.^2,[3 4])) ;
        surf(xp(:,:,1),xp(:,:,2),0*u(:,:,1),clr,'edgecolor','none')
        colorbar
    case 3
        vol = volshow(kappa) ;
        vol.RenderingStyle = 'CinematicRendering' ;
        viewer = vol.Parent ;
        viewer.BackgroundGradient = 'off' ;
        viewer.BackgroundColor = 'w' ;
end


return
%% TESTS
clc
disp('TESTS')
disp('Fourier Transforms')
u = randn([N D])+1i*randn([N D]) ;
U = to_fourier(u,N) ;
uu = to_spatial(U,N) ; 
norm(uu(:)-u(:))
disp('Tensprod-Energy')
eps = randn([N D D]) ;
sig = randn([N D D]) ;
tensprod(sig,eps,2) ; 
size(ans)-N
disp('Tensprod-stresses')
C = randn([N D D D D]) ;
sig = tensprod(C,eps,2) ;
size(sig)-[N D D]
disp('Tensprod-uniform stresses') ;
C0 = mean(C,1:D) ;
sig = tensprod(C,eps,2) ;
size(sig)-[N D D]
disp('Derivation-grad(u)')
u = randn([N D]) ;
F = grad(u,xi) ;
size(F)-[N D D]
disp('Derivation-div(u)')
u = randn([N D]) ;
d = div(u,xi) ;
size(d)-[N]
disp('Derivation-div(sig)')
sig = randn([N D D]) ;
d = div(sig,xi) ;
size(d)-[N D]
disp('assemble_isoC-full')
kappa = randn([N 1]) ;
mu = randn([N 1]) ;
tic ; C = assemble_isoC(kappa,mu) ; toc
size(C)-[N D D D D]
disp('C*eps')
eps = randn([N D D]) ;
tic ; sig = tensprod(C,eps,2) ; toc
size(sig)-[N D D]
disp('apply_isoC-full')
tic ; sig_app = apply_isoC(eps,kappa,mu) ; toc
norm(sig(:)-sig_app(:))
disp('P0_assembly')
NABLA = assemble_nabla(xi) ;
tic ; P0 = tensprod(conj(permute(NABLA,[1:D D+3 D+1 D+2])),tensprod(C0,NABLA,2),2) ; toc
tic ; P0_direct = div(tensprod(C0,grad(reshape(eye(D),[ones(1,D) D D]),xi),2),xi) ; toc
norm(P0(:)-P0_direct(:))


disp('gradient of known function')
k = 2*pi./L.*ones(1,D) ;
f = exp(1i*sum(reshape(k,[ones(1,D) D]).*x,D+1)) ;
gf = to_spatial(grad(to_fourier(f,N),xi),N) ;
(gf - 1i*reshape(k,[ones(1,D) D]).*f) ; norm(ans(:))
lf = to_spatial(div(to_fourier(gf,N),xi),N) ;
(lf - sum(k.^2).*f) ; norm(ans(:))

NABLA = assemble_nabla(xi,false,'willot') ;


%% FUNCTIONS

function x = grid2vec(x,N)
    x = reshape(x,[prod(N) size(x,numel(N)+1:ndims(x))]) ;
end

function x = vec2grid(x,N)
    x = reshape(x,[N size(x,numel(N)+1:ndims(x))]) ;
end

function x = grid2dof(x)
    x = x(:) ;
end

function x = dof2grid(x,N,n)
    if nargin<3 
        n = numel(x)/prod(N) ; 
        D = numel(N) ;
        if mod(log(n)/log(D),1)<sqrt(eps)
            n = D*ones(1,D) ;
        end
    end
    x = reshape(x,[N n]) ;
end

function x = to_fourier(x,N)
    for d = 1:numel(N)
        x = fft(x,N(d),d) ;
        % x = fftshift(x,d) ;
    end
end

function x = to_spatial(x,N)
    for d = 1:numel(N)
        % x = ifftshift(x,d) ;
        x = ifft(x,N(d),d) ;
    end
end

function x = shift_fft(x,N,inv)
    if nargin<3 ; inv = false ; end
    for d = 1:numel(N)
        if inv
            x = ifftshift(x,d) ;
        else
            x = fftshift(x,d) ;
        end
        
    end
end

function C = tensprod(A,B,n,D)
    if nargin<3 ; n=0 ; end
    if nargin<4 ; D = max(size(A,ndims(A)),size(B,ndims(B))) ; end
    N = max(size(A,1:D),size(B,1:D)) ;
    nA = ndims(A)-D ;
    C = A.*reshape(B,[size(B,1:D) ones(1,nA-n) size(B,(D+1):ndims(B))]) ;
    if n>0
        C = sum(C,ndims(A)-n+(1:n)) ;
    end
    C = reshape(C,[N size(A,D+(1:nA-n)) size(B,D+n+1:ndims(B))]) ;
end

function dx = apply_nabla(x,xi,conjugate,method)
    if nargin<3 ; conjugate = false ; end
    if nargin<4 ; method = 'analytic' ; end
    switch method
        case 'analytic'
            der = 1i*xi ;
            if conjugate ; dx = tensprod(conj(der),x,1) ;
            else ; dx = tensprod(der,x,0) ;
            end
        case 'willot'
            D = size(xi,ndims(xi)) ; N = size(xi,1:D) ;
            switch D
                case 2
                    dxi = xi(2,2,:)-xi(1,1,:) ;
                    dx = 2*pi./N./dxi ;
                    der = cat(D+1 ...
                                ,(2i/dx(1))*sin(.5*dx(1).*xi(:,:,1)).*cos(.5*dx(2).*xi(:,:,2)) ...
                                ,(2i/dx(2))*cos(.5*dx(1).*xi(:,:,1)).*sin(.5*dx(2).*xi(:,:,2)) ...
                            ) ;
                case 3
                    error
            end
            if conjugate ; dx = tensprod(conj(der),x,1) ;
            else ; dx = tensprod(der,x,0) ;
            end
        otherwise
            error("Unknowm derivation method")
    end
end

function NABLA = assemble_nabla(xi,varargin)
    D = ndims(xi)-1 ;
    x = reshape(eye(D),[ones(1,D) D D]) ;
    NABLA = apply_nabla(x,xi,varargin{:}) ;
end

function g = grad(x,xi,method)
    if nargin<3 ; method = 'analytic' ; end
    g = apply_nabla(x,xi,false,method) ;
end

function d = div(x,xi,method)
    if nargin<3 ; method = 'analytic' ; end
    d = apply_nabla(x,xi,true,method) ;
end

function eps = sym(g)
    D = size(g,ndims(g)) ;
    eps = reshape(g,[],D,D) ;
    eps = .5*(eps+permute(eps,[1 3 2])) ;
    eps = reshape(eps,[size(g,1:ndims(g)-2) D D]) ;
end

function C = assemble_isoC(kappa,mu,D)
    if nargin<3 ; D = ndims(kappa+mu) ; end
    ii = reshape(1:D,[ones(1,D) D 1 1 1]) ;
    jj = reshape(1:D,[ones(1,D) 1 D 1 1]) ;
    kk = reshape(1:D,[ones(1,D) 1 1 D 1]) ;
    ll = reshape(1:D,[ones(1,D) 1 1 1 D]) ;
    I = (1/2)*( (ii==kk).*(jj==ll) + (ii==ll).*(jj==kk) ) ; % identity: I:e = e
    J = (1/D)*(ii==jj).*(kk==ll) ; % spherical part: J:e = (1/D)*Tr(e)*eye(D)
    K = I-J ; % deviatoric part
    C = 3*kappa.*J + 2*mu.*K ;
end

function sig = apply_isoC(eps,kappa,mu) 
    D = size(eps,ndims(eps)) ;
    N = max(size(eps,1:D),size(kappa+mu,1:D)) ;
    eps = reshape(eps,[],D,D) ;
    % trace
    t = eps(:,1,1) ;
    for d = 2:D 
        t = t + eps(:,d,d) ;
    end
    % spherical part
    sph = (1/D)*t.*reshape(eye(D),[1 D D]) ;
    % deviatoric part
    dev = .5*(eps + permute(eps,[1 3 2])) - sph ;
    % stresses
    sig = 3*kappa(:).*sph + 2*mu(:).*dev ;
    sig = reshape(sig,[N D D]) ;
end










