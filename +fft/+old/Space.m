classdef Space < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
% A tensor space defined on a regular grid
    

    properties
        Grid pkg.fft.Grid
    end

    properties (Dependent)
        Order(1,1) {mustBeInteger, mustBeNonnegative}
    end

    properties 
        Size(1,:) {mustBeInteger, mustBeNonnegative} % Size=[] is a scalar field
        Fourier(1,1) logical = false ; % is the field in fourier domain ?
    end

    properties (Hidden,SetAccess=protected)
        DataTemplate = 1 % used to keep trace of the data type in pkg.fft.Operators
    end
    
    methods
        function this = Space(varargin)
        % Constructor
        % this = pkg.fft.Space(space) ; % transparent
        % this = pkg.fft.Space(variable) ; % extract the space from a pkg.fft.Variable
        % this = pkg.fft.Space(grid) ; % scalar space on a given pkg.fft.Grid
        % this = pkg.fft.Space(...,Name,Value) ;
            if nargin==0 ; return ; end
            switch class(varargin{1})
                case 'pkg.fft.Space'
                    for prop = properties(pkg.fft.Space())'
                        this.(prop{:}) = varargin{1}.(prop{:}) ;
                    end
                case 'pkg.fft.Variable'
                    for prop = properties(pkg.fft.Space())'
                        this.(prop{:}) = varargin{1}.(prop{:}) ;
                    end
                case 'pkg.fft.Grid'
                    this.Grid = varargin{1} ;
                otherwise % guess a "constant" space from some numeric data
                    this.Grid = pkg.fft.Grid([]) ; % scalar grid
                    this.Size = size(varargin{1}) ;
                    this.DataTemplate = varargin{1}(1) ;
            end
            for arg = 2:2:numel(varargin)
                this.(varargin{arg}) = varargin{arg+1} ;
            end
        end

        function set.Size(this,sz)
            this.Size = sz(1:find(sz>1,1,'last')) ; % remove trailing singleton dimensions
        end

        function set.Order(this,order)
            this.Size = [] ;
            if order~=0 % not scalar
                this.Size = repmat(this.Grid.ndims,[1 order]) ;
            end
        end

        function order = get.Order(this)
            order = numel(this.Size) ;
        end

        function sz = datasize(this) 
            sz = [this.Grid.N this.Size] ;
        end

        function n = datalength(this) 
            n = prod(this.datasize) ;
        end

        function field = setdata(this,data)
            field = pkg.fft.Variable(this,'Data',data) ;
        end
    
    end


    methods

        function x = datacast(x,fcn)
            x = copy(x) ;
            x.DataTemplate = fcn(x.DataTemplate) ;
        end
        function x = single(x) ; x = datacast(x,@single) ; end
        function x = double(x) ; x = datacast(x,@double) ; end
        function x = complex(x) ; x = datacast(x,@complex) ; end
        function x = gpuArray(x) ; x = datacast(x,@gpuArray) ; end

        function x = fouriertransform(x,inverse)
        % (Inverse) Fourier Transformed space
            if nargin<2 ; inverse = false ; end % do normal FFT by default
            if x.Fourier ~= inverse ; return ; end % the field is already in the queried space
            x = copy(x) ;
            x.Fourier = ~inverse ;
            x.DataTemplate = complex(x.DataTemplate) ;
        end
        function x = fft(x) ; x = fouriertransform(x,false) ; end
        function x = ifft(x) ; x = fouriertransform(x,true) ; end

        function C = arrayfun(fun,varargin)
        % convert all inputs to spaces
            isspace = cellfun(@(arg)isa(arg,'pkg.fft.Space'),varargin) ;
            if ~all(isspace)
                varargin = cellfun(@pkg.fft.Space,varargin) ;
            else
                varargin = [varargin{:}] ;
            end
            if range([varargin(isspace).Fourier])~=0 ; warning('Space arguments with different Fourier domains !') ; end
        % build the common grid
            grid = varargin(1).Grid ;
            for arg = 2:numel(varargin) 
                grid = grid.*varargin(arg).Grid ; 
            end
            C = pkg.fft.Space(grid) ;
        % guess the data sizes
            order = max([varargin.Order]) ;
            sz = cellfun(@(sz)[sz ones(1,order-numel(sz))],{varargin.Size},'uni',false) ;
            C.Size = max(cat(1,sz{:}),[],1) ;
            C.Fourier = any([varargin(isspace).Fourier]) ;
        % Data type template
            C.DataTemplate = fun(varargin.DataTemplate) ;
        end
        function varargout = plus(varargin) ; [varargout{1:nargout}] = arrayfun(@plus,varargin{:}) ; end
        function varargout = uplus(varargin) ; [varargout{1:nargout}] = arrayfun(@uplus,varargin{:}) ; end
        function varargout = minus(varargin) ; [varargout{1:nargout}] = arrayfun(@minus,varargin{:}) ; end
        function varargout = uminus(varargin) ; [varargout{1:nargout}] = arrayfun(@uminus,varargin{:}) ; end
        function varargout = times(varargin) ; [varargout{1:nargout}] = arrayfun(@times,varargin{:}) ; end
        function varargout = real(varargin) ; [varargout{1:nargout}] = arrayfun(@real,varargin{:}) ; end
        function varargout = imag(varargin) ; [varargout{1:nargout}] = arrayfun(@imag,varargin{:}) ; end
        function varargout = abs(varargin) ; [varargout{1:nargout}] = arrayfun(@abs,varargin{:}) ; end
        function varargout = conj(varargin) ; [varargout{1:nargout}] = arrayfun(@conj,varargin{:}) ; end
        function varargout = sign(varargin) ; [varargout{1:nargout}] = arrayfun(@sign,varargin{:}) ; end
        function varargout = power(varargin) ; [varargout{1:nargout}] = arrayfun(@power,varargin{:}) ; end

        function x = transpose(x)
            x = copy(x) ;
            switch x.Order
                case 0 % do nothing
                case 1 % column vector-> row matrix
                    x.Size = [1 x.Size] ;
                otherwise % second order and more-> swap the two first dimensions
                    perm = [2 1 3:x.Order] ;
                    x.Size = x.Size(perm) ;
            end
        end

        function x = ctranspose(x)
            x = conj(transpose(x)) ;
        end

        function C = tensorprod(A,B,order)
        % TENSORPROD of two gridded tensor field variables at a given order
            if isa(A,'pkg.fft.Space') && isa(B,'pkg.fft.Space') && A.Fourier~=B.Fourier
                warning('Mixing spaces with different Fourier domains !') ;
            end
            A = pkg.fft.Variable(A) ; B = pkg.fft.Variable(B) ;
            C = pkg.fft.Variable(A.Grid.*B.Grid,'Size',[A.Size(1:end-order) B.Size(order+1:end)]) ;
            C.Fourier = A.Fourier || B.Fourier ;
            C.DataTemplate = A.DataTemplate*B.DataTemplate ;
        end

        function g = grad(x,nabla)
            if nargin<2 ; nabla = pkg.fft.fields.nabla(x.Grid) ; end
            g = tensorprod(nabla,fft(x),0) ;
        end
            
        function d = div(x,nabla)
            if nargin<2 ; nabla = pkg.fft.fields.nabla(x.Grid) ; end
            d = tensorprod(conj(nabla),fft(x),1) ;
        end

        function x = sym(x)
            x = 0.5.*(x+x.') ;
        end

        function t = trace(x)
            if x.Order<2 ; error('TRACE applies to tensors of order>2 !') ; end
            if x.Size(1)~=x.Size(2) ; error('TRACE applies to tensors having the two first dimensions with the same size !') ; end 
            t = pkg.fft.Space(x.Grid,'Size',x.Size(3:end),'Fourier',x.Fourier) ;
        end

        function x = gridfun(x,~,dims)
            x = copy(x) ;
            if nargin<3 ; dims = 1:x.Grid.ndims ; end
            x.Grid = contract(x.Grid,dims) ;
        end
        function varargout = gridmean(x,varargin) ; [varargout{1:nargout}] = gridfun(x,@mean,varargin{:}) ; end
        function varargout = gridmedian(x,varargin) ; [varargout{1:nargout}] = gridfun(x,@median,varargin{:}) ; end
        function varargout = gridsum(x,varargin) ; [varargout{1:nargout}] = gridfun(x,@sum,varargin{:}) ; end

        function x = tensorfun(x,~,dims)
            x = copy(x) ;
            if nargin<3 ; dims = 1:x.Order ; end
            x.Size(dims) = 1 ; % this might clear the data
        end
        function varargout = tensormean(x,varargin) ; [varargout{1:nargout}] = tensorfun(x,@mean,varargin{:}) ; end
        function varargout = tensormedian(x,varargin) ; [varargout{1:nargout}] = tensorfun(x,@median,varargin{:}) ; end
        function varargout = tensorsum(x,varargin) ; [varargout{1:nargout}] = tensorfun(x,@sum,varargin{:}) ; end
        function varargout = tensormin(x,y,varargin) 
            if nargin<2 || isempty(y)
                if x.Order==0
                    [varargout{1:nargout}] = copy(x,nargout) ; 
                    return ;
                else
                    y = [] ; 
                end 
            end
            [varargout{1:nargout}] = tensorfun(x,@(x,dims)min(x,y,dims),varargin{:}) ; 
        end

        function x = tensornorm(x)
            x = copy(x) ;
            x.Order = 0 ;
        end

        function x = tensorinv(x)
            if x.Order~=2 ; error('INV applies to second-order tensors only !') ; end
        end


    end

end

