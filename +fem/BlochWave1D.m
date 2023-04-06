%% TEST 1D BLOCH-WAVE PROPAGATION, HELMHOLTZ MEDIA
% A.d²u/dx² = d²u/dt²
% u(x,y) = U(x).exp(1i*k*y)
clc,clearvars

L = 1 ; % Unit cell length
nX = 100 ; % number of elements/unit cell
Fa = [1 .1 .0*exp(2i*pi/180*0) .1 0 0 0 0 0 0] ; % fourier expansion of impedance A
nK = 100 ; % wavenumber discretization
nModes = min(10,nX) ;


x = linspace(0,L,nX+1)' ;
mesh = pkg.geometry.mesh.Mesh(x) ;

A = ifft([Fa zeros(1,nX-2*numel(Fa)+1) flip(conj(Fa(2:end)))])*nX ;
clf ; area(mesh.centroid,A) ; axis tight

T = [speye(nX) ; speye(1,nX)] ;

[ee,we,ie] = mesh.integration() ;
N = mesh.interpMat(ee,ie)*T ;
D = mesh.gradMat(1,ee,ie)*T ;

Aw = diag(sparse(we(:).*A(:))) ;
K0 = D'*Aw*D ;
K1 = D'*Aw*N - N'*Aw*D;
K2 = N'*Aw*N ;
M = N'*diag(sparse(we))*N ;

k = linspace(-pi,pi,nK)/L ;
w = zeros(min(nModes,nX),nK) ;
U = zeros([nX size(w)]) ;
for ii = 1:nK
    K = K0 + 1i*K1*k(ii) + K2*k(ii)^2 ;
    if nModes>=nX
        [Ui,w2] = eig(full(K),full(M),'vector') ;
        [~,is] = sort(real(w2),'ascend') ;
        U(:,:,ii) = Ui(:,is(1:nModes)) ;
        w(:,ii) = sqrt(w2(is(1:nModes))) ;
    else
        [U(:,:,ii),w2] = eigs(K,M,nModes,'sm') ;
        w(:,ii) = sqrt(diag(w2)) ;
    end
end

c = w./k ; 
wc = pi*mean(sqrt(A))/L ;
w = w/wc ; 

clf ; plot(real(w),real(k),'.k')



