clearvars

% CPU/GPU computations
datacast = @(data)gpuArray(single(data)) ;
% Unit Cell dimensions
L = [1 1 1] ;
% discrete grid
N = (2^7)*ones(size(L)) ;
grid = datacast(pkg.fft.old.Grid(N,'dx',L./N)) ;
D = double(ndims(grid)) ;

% Material phase field
contrast = 1e3 ;
x0 = .5*L ; R0 = .3*min(L) ; % soft pores
xS = [0*L;L;eye(D).*L;(1-eye(D)).*L] ; Rs = R0 ; % rigid inclusions
x = grid.coordinates() ;
r20 = tensorsum((x-x0').^2,1) ;
phi0 = tensormin(r20-(R0.^2)',[]) ;
r2s = tensorsum((x-xS').^2,1) ;
phiS = tensormin(r2s-(Rs.^2)',[]) ;
phi = .5.*(sign(phi0)-sign(phiS)) ;

% material properties
kappa0 = 1 ; mu0 = 1 ; rho0 = 1 ;
modulation = sqrt(contrast).^phi ;
kappa = kappa0.*modulation ;
mu = mu0.*modulation ;
rho = rho0.*modulation ;
if D==2 ; clf ; axis equal tight ; plot(modulation) ; set(gca,'colorscale','log') ; colorbar ; end

% Displacement field
u = pkg.fft.old.Field(grid,'Order',1,'Fourier',true) ; 

% Derivation operator
method = 'willot' ;
nabla = pkg.fft.old.fields.nabla(grid,method) ;
% EPS = @(u)ifft(sym(grad(u,nabla))) ;
EPS = @(u)ifft(symgrad(u,nabla)) ;
Id = datacast(pkg.fft.old.Field((1/D)*eye(D))) ;
sph = @(EPS)tensorprod(Id,trace(EPS),0) ;
SIG = @(kappa,mu,EPS)((3.*kappa)-(2.*mu)).*sph(EPS) + (2.*mu).*EPS ;
f = @(SIG)fft(div(SIG,nabla)) ;

E = .1*[0 1 ; 1 0] ; E = padarray(E,[D-2 D-2],0,'post') ;
b = -div(SIG(kappa,mu,E),nabla) ;
if D==2 ; clf ; axis equal tight ; pl = plot(real(ifft(b))) ; colorbar ; end

A0 = assemble(@(u)f(SIG(kappa0,mu0,EPS(u))),u) ;
singular = all(abs(A0.Data)<sqrt(eps),1:2) ;
A0.Data = A0.Data + singular.*eye(D) ;
iA0 = inv(A0) ;
iA0_op = @(x)fft(tensorprod(iA0,u.setdata(x),1)).getdata(true) ;

% C = assemble(@(EPS)SIG(kappa,mu,EPS),EPS(u)) ;
% SIG = @(kappa,mu,EPS)tensorprod(C,EPS,2) ;

A = @(u)f(SIG(kappa,mu,EPS(u))) ;
A_op = @(x)fft(A(u.setdata(x))).getdata(true) ;

profile on
tic ; u.Data = pcg(A_op,b.getdata(true),1e-3,100,iA0_op) ; toc
profile off
u = ifft(u) ;

x = grid.coordinates() ;
e = ifft(EPS(u)) + E ;
u = u + tensorprod(pkg.fft.old.Field(E),x,1) ;

if D==2 ; clf ; axis equal tight off ; pl = plot(x+real(u),norm(e)) ; colorbar ; end

return

%% Timing...

%profile on
disp('=== Timing ===')

% Tensorprod function
disp('Tensor product')
EPS = pkg.fft.old.Field(grid,'Order',2) ;
C = pkg.fft.old.Field(grid,'Order',4) ;
fcn = @()tensorprod(C,EPS,2) ;
CPU_GPU_compare(fcn,[C EPS]) ;

% Hooke's law
disp('Hooke''s law')
EPS = pkg.fft.old.Field(grid,'Order',2) ;
kappa = pkg.fft.old.Field(grid,'Order',0) ;
mu = copy(kappa) ;
fcn = @()SIG(kappa,mu,EPS) ;
CPU_GPU_compare(fcn,[kappa mu EPS]) ;

%% FFT/IFFT
disp('FFT/IFFT')
u = pkg.fft.old.Field(grid,'Order',1) ;
fcn = @()ifft(fft(u)) ;
CPU_GPU_compare(fcn,u) ;


% profile off
% profile viewer


%% UTILS
function A = assemble(A_op,input_field)
    sz = input_field.Size ;
    Id = pkg.fft.old.Field(input_field.Grid,'Size',[sz sz],'Fourier',input_field.Fourier) ;
    Id.Data = repmat(eye(prod(sz)),[1 1 Id.Grid.N]) ;
    A = A_op(Id) ;
end


function [CPUtime,GPUtime] = CPU_GPU_compare(fcn,fields,dataclass)
    if nargin<3 ; dataclass = 'single' ; end
    % fields = arrayfun(@copy,fields) ;
    fields = arrayfun(@(field)field.setdata(@(sz)randn(sz,dataclass)+1i*randn(sz,dataclass)),fields) ;
    % Push data into MAIN memory
    fields = arrayfun(@(field)field.setdata(gather(field.getdata)),fields) ; 
    % CPU version
    CPUtime = timeit(fcn) ; 
    disp("  CPU:"+string(CPUtime))
    % Push data to GPU
    fields = arrayfun(@(field)field.setdata(gpuArray(field.getdata)),fields) ; 
    % check that the output is still on the GPU
    if ~strcmp(class(fcn().Data),'gpuArray') ; warning('output not on GPU') ; end
    % GPU version
    GPUtime = gputimeit(fcn) ; 
    disp("  GPU:"+string(GPUtime))
end

