function var = nabla(var,method,conjugate)
%NABLA derivation of a gridded tensor field variable

    if nargin<2 || isempty(method) ; method = 'analytic' ; end
    if nargin<3 ; conjugate = false ; end

% Grid step & Fourier Coordinates
    if var.Grid.Fourier
        xi = var.Grid.coordinates() ;
        dx = ifft(var.Grid).dx ;
    else
        dx = var.Grid.dx ;
        var = pkg.fft.old.functions.fft(var) ;
        xi = var.Grid.coordinates() ;
    end
% Derivation kernel
    kernel = copy(xi) ;
    switch method
        case 'analytic'
            kernel.Data = 1i*xi.Data ;
        case 'willot'
            kernel.Data(:,:) = .5i*cast(dx,class(kernel.Data))*kernel.Data(:,:) ;
            idx = cast(dx\eye(size(dx,1)),class(kernel.Data)) ;
            kernel.Data(:,:) = idx*(exp(kernel.Data(:,:))-2*(1-2*eye(size(dx,1)))*exp(-kernel.Data(:,:))) ;
    end
    if conjugate ; kernel.Data = conj(kernel.Data) ; end
% Apply the kernel
    if conjugate % divergence of the tensor field
        var = pkg.fft.old.functions.tensorprod(var,kernel,1) ;
    else % gradient of the tensor field
        var = pkg.fft.old.functions.tensorprod(var,kernel,0) ;
    end


end

