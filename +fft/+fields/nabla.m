function kernel = nabla(grid,method)
% Fourier-domain derivation operator on a gridded tensor field

    if nargin<2 || isempty(method) ; method = 'analytic' ; end

% Grid step & Fourier Coordinates
    xi = grid.fouriercoordinates ;
    dx = grid.dx ;

% Derivation FFT-kernel
    kernel = pkg.fft.Field(grid,'Order',1,'Fourier',true) ;
    switch method
        case 'analytic' % fft(df_dx) = 1i*xi*fft(f)
            kernel.Data = 1i*xi.Data ;
        otherwise % discrete formulations
            xi = reshape(xi.Data,[],xi.Size) ;
            xi_dx = xi*dx.' ;
            inv_dx = dx\eye(size(dx,1),class(xi)) ;
            switch method
                case 'finitedifferences' % df_dx = (f(x+dx/2)-f(x-dx/2))/dx
                    ker = sin(.5*xi_dx)*(2i*inv_dx.') ;
                case 'willot' % mean derivative over the neightbor nodes
                    ker = sin(.5*xi_dx)*(2i*inv_dx.') ;
                    for d = 1:kernel.Grid.ndims
                        od = [1:d-1 d+1:kernel.Grid.ndims] ;
                        ker(:,od) = ker(:,od).*cos(.5*xi_dx(:,od)) ;
                    end
            end
            kernel.Data = ker ;
    end

% Remove "ambiguous" wavevectors where xi=\pm N*pi/L
    for d = find(mod(kernel.Grid.N,2)==0) % for any dimension with even number of fourier frequencies
        selvec = 1-full(sparse(1,double(kernel.Grid.N(d))/2+1,1,1,kernel.Grid.N(d))) ; 
        selvec = reshape(selvec,[ones(1,d-1) kernel.Grid.N(d) 1 1]) ;
        kernel.Data = kernel.Data.*selvec ;
    end

end

