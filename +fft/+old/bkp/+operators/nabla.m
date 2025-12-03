classdef nabla < pkg.fft.old.operator
%NABLA partial derivative operator
    
    properties
        Method char = 'analytic' ;
        Conjugate(1,1) logical = false ;
        Operator(1,1) pkg.fft.old.operator
    end
    
    methods
        function this = nabla(f,method,conjugate)
        % Constructor
            this.Inputs = f ;
            if nargin>1 && ~isempty(method) ; this.Method = method ; end
            if nargin>2 ; this.Conjugate = conjugate ; end
        % Grid step & Fourier Coordinates
            if this.Inputs.fourier
                xi = this.Inputs.Grid.coordinates() ;
                dx = ifft(this.Inputs.Grid).dx ;
            else
                dx = this.Inputs.Grid.dx ;
                xi = fft(this.Inputs.Grid).coordinates() ;
            end
        % Derivation kernel
            kernel = copy(xi) ;
            switch this.Method
                case 'analytic'
                    kernel.Data = 1i*xi.Data ;
                case 'willot'
                    kernel.Data(:,:) = .5i*cast(dx,class(kernel.Data))*kernel.Data(:,:) ;
                    idx = cast(dx\eye(size(dx,1)),class(kernel.Data)) ;
                    kernel.Data(:,:) = idx*(exp(kernel.Data(:,:))-2*(1-2*eye(size(dx,1)))*exp(-kernel.Data(:,:))) ;
            end
            if this.Conjugate ; kernel.Data = conj(kernel.Data) ; end
        % Operator
            if this.Conjugate % divergence of the tensor field
                this.Operator = pkg.fft.old.operators.tensorprod(kernel,this.Inputs,1) ;
            else % gradient of the tensor field
                this.Operator = pkg.fft.old.operators.tensorprod(kernel,this.Inputs,0) ;
            end
        % Output field
            this.Outputs = this.Operator.Outputs ;
        end

        function apply(this)
            this.Operator.apply() ;
        end


    end
end

