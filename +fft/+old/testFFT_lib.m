clc
clearvars

L = single([1 1]) ;
N = 128 ; %*ones(size(L)) ;

% Grid creation
grid = pkg.fft.old.Grid(N,'dx',L./N) ;
% Grid composition
xgrid = pkg.fft.old.Grid([N 1],'dx',L(1:2)./N) ;
ygrid = pkg.fft.old.Grid([1 N],'dx',L(1:2)./N) ;
tp_grid = xgrid.*ygrid ;

% Coordinates on grids
x = grid.coordinates() ;
xi = grid.fouriercoordinates() ;

% Fields of different orders
% Pressure field
p = pkg.fft.old.Field(grid,'Order',0) ;
% Displacement field
u = pkg.fft.old.Field(grid,'Order',1) ;
% Strain field
EPS = pkg.fft.old.Field(grid,'Order',2) ;
% Stiffness tensor field
C = pkg.fft.old.Field(grid,'Order',4) ;
% [displacement+rotation] field
TORS = pkg.fft.old.Field(grid,'Size',[3 2]) ;

% Fill data
u.Data = @(n)randn(n)+1i*randn(n) ;
C.Data = @randn ;

% Functions
% Gradient
EPS = sym(grad(u)) ; % /!\ in Fourier domain !
% Tensor product
SIG = tensorprod(C,ifft(EPS),2) ;
% Any other function
kappa = pkg.fft.old.Field(grid).setdata(@rand) ;
mu = pkg.fft.old.Field(grid).setdata(@rand) ;
sph_EPS = ((1/double(grid.ndims))*eye(grid.ndims)).*trace(EPS) ;
dev_EPS = EPS-sph_EPS ;
SIG_iso = (3.*kappa).*ifft(sph_EPS) + (2.*mu).*ifft(dev_EPS) ;
% Divergence
f = div(SIG) ;

return

%% Operator Assembly
method = 'analytic' ;
nabla = pkg.fft.old.fields.nabla(grid,method) ;
C0 = gridmean(C,1:C.Grid.ndims) ;
A_op = @(u)div(tensorprod(C,sym(grad(u,nabla)),2),nabla) ;
Id = pkg.fft.old.Field(grid,'Order',2) ;
Id.Data = repmat(eye(Id.Grid.ndims),[1 1 Id.Grid.N]) ;
tic ; A = A_op(Id) ; toc

return



%% TESTS
% FFT computations
u.Data = complex(single(randn([u.Size u.Grid.N]))) ;
disp('---Direct CPU FFT')
tic ; Fu = fft(u) ; toc
disp("|Fu-u|="+string(norm(u.Data(:)-Fu.Data(:))))
disp('---Direct GPU FFT')
u.Data = gpuArray(u.Data) ;
tic ; Fu = fft(u) ; toc
disp("|Fu-u|="+string(norm(u.Data(:)-Fu.Data(:))))
disp('---In-Place GPU FFT')
u.Data = gpuArray(u.Data) ;
tic ; fft(u) ; toc
disp("|Fu-u|="+string(norm(u.Data(:)-Fu.Data(:))))
%% function derivation
disp('---gradient of known function')
f = pkg.fft.old.field(grid,'Order',0) ;
k = 2*pi./L ;
f.Data = reshape(exp(1i*k*x.Data(:,:)),f.Grid.N) ;
gradF = pkg.fft.old.operators.nabla(fft(f)) ;
gradF.apply() ;
laplacianF = pkg.fft.old.operators.nabla(gradF.Outputs,[],true) ;
laplacianF.apply() ;
Lf = ifft(laplacianF.Outputs) ;
%%
gf = to_spatial(grad(to_fourier(f,N),xi),N) ;
(gf - 1i*reshape(k,[ones(1,D) D]).*f) ; norm(ans(:))
lf = to_spatial(div(to_fourier(gf,N),xi),N) ;
(lf - sum(k.^2).*f) ; norm(ans(:))





