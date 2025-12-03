classdef Field < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
% A tensor field defined on a regular grid
    

    properties
        Grid pkg.fft.Grid
    end

    properties (Dependent)
        Order(1,1) {mustBeInteger, mustBeNonnegative}
    end

    properties 
        Size(1,:) {mustBeInteger, mustBeNonnegative} % Size=[] is a scalar field
        Fourier(1,1) logical = false ; % is the field data in fourier domain ?
        Data
    end
    
    methods
        function this = Field(varargin)
        % Constructor
        % this = field(field) ; % transparent
        % this = field(data) ; % conversion to a constant field
        % this = field(grid) ;
        % this = field(...,Name,Value) ;
            if nargin==0 ; return ; end
            switch class(varargin{1})
                case 'pkg.fft.Field'
                    this = varargin{1} ;
                case 'pkg.fft.Grid'
                    this.Grid = varargin{1} ;
                otherwise
                    this.Grid = pkg.fft.Grid([]) ; % scalar grid
                    this.Size = size(varargin{1}) ;
                    this.Data = varargin{1} ;
            end
            for arg = 2:2:numel(varargin)
                this.(varargin{arg}) = varargin{arg+1} ;
            end
        end

        function set.Size(this,sz)
            this.Size = sz(1:find(sz>1,1,'last')) ; % remove trailing singleton dimensions
            this.Data = [] ;
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

        function data = get.Data(this)
            data = this.Data ;
        end

        function set.Data(this,data)
            if isempty(data) ; this.Data = [] ; return ; end
            switch class(data)
                case 'function_handle'
                    data = data(this.datasize) ;
                otherwise
            end
            this.Data = reshape(data,[this.datasize 1 1]) ;
        end

        function this = setdata(this,data)
            this.Data = data ;
        end

        function data = getdata(this,flatten)
            data = this.Data ;
            if nargin>1 && flatten ; data = data(:) ; end
        end

        function this = softcopy(this,do_copy) % might do in-place computations if possible
            if nargin<2 || do_copy ; this = copy(this) ; end
        end

        function pl = plot(varargin)
            if nargin>1 && isa(varargin{2},'pkg.fft.Field') % plot(x,field,...)
                x = varargin{1} ;
                field = varargin{2} ;
                varargin = varargin(3:end) ;
            else % plot(field,...)
                field = varargin{1} ;
                x = field.Grid.coordinates ;
                varargin = varargin(2:end) ;
            end
            if field.Fourier
                field = copy(field) ; % to not modify the plotted field data
                x = field.Grid.fouriercoordinates ;
                for d = 1:field.Grid.ndims
                    field.Data = fftshift(field.Data,d) ;
                    x.Data = fftshift(x.Data,d) ;
                end
            end
            pl = plot(field.Grid,varargin{:}) ;
            pl.FaceColor = 'interp' ;
            pl.EdgeColor = 'none' ;
            pl.XData = gather(x.Data(:,:,1)) ;
            pl.YData = gather(x.Data(:,:,2)) ;
            if field.Order==0
                pl.CData = gather(field.Data) ;
            else
                pl.CData = gather(norm(field).Data) ;
            end
        end
    
    end


    methods

        function x = single(x)
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            x.Data = single(x.Data) ;
        end

        function x = double(x)
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            x.Data = double(x.Data) ;
        end

        function x = complex(x)
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            x.Data = complex(x.Data) ;
        end

        function x = gpuArray(x)
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            x.Data = gpuArray(x.Data) ;
        end

        function x = fouriertransform(x,inverse)
            if nargin<2 ; inverse = false ; end % do FFT by default
            if x.Fourier~=inverse ; return ; end % the field is already in the queried space
            if ~isempty(x.Data)
                for d = 1:x.Grid.ndims
                    if inverse
                        x.Data = ifft(x.Data,x.Grid.N(d),d) ;
                    else
                        x.Data = fft(x.Data,x.Grid.N(d),d) ;
                    end
                end
            end
            x.Fourier = ~inverse ;
        end

        function x = fft(x)
        % FFT of a gridded tensor field
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            x = fouriertransform(x,false) ;
        end

        function x = ifft(x)
        % inverse FFT of a gridded tensor field
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            x = fouriertransform(x,true) ;
        end

        function C = elemfun(fun,varargin)
        % convert all inputs to fields
            isfield = cellfun(@(arg)isa(arg,'pkg.fft.Field'),varargin) ;
            if ~all(isfield)
                varargin = cellfun(@pkg.fft.Field,varargin) ;
            else
                varargin = [varargin{:}] ;
            end
            if range([varargin(isfield).Fourier])~=0 ; warning('Field arguments with different Fourier domains !') ; end
        % build the common grid
            grid = varargin(1).Grid ;
            for arg = 2:numel(varargin) 
                grid = grid.*varargin(arg).Grid ; 
            end
            C = pkg.fft.Field(grid) ;
        % match the data sizes
            order = max([varargin.Order]) ;
            if numel(unique([varargin.Grid]))==1 % all fields are on the same grid
                data = fun(varargin.Data) ; 
            else
                data = arrayfun(@(field)reshape(field.Data,[field.Grid.N ones(1,ndims(C.Grid)-ndims(field.Grid)) field.Size]),varargin,'uni',false) ;
                data = fun(data{:}) ;
            end
            C.Size = size(data,ndims(C.Grid)+(1:order)) ;
            C.Data = data ;
            C.Fourier = any([varargin(isfield).Fourier]) ;
        end

        function C = tensorprod(A,B,order)
        % TENSORPROD of two gridded tensor field variables at a given order
            if isa(A,'pkg.fft.Field') && isa(B,'pkg.fft.Field') && A.Fourier~=B.Fourier
                warning('Mixing two grids that are not on different spaces !') ;
            end
            A = pkg.fft.Field(A) ; B = pkg.fft.Field(B) ;
            C = pkg.fft.Field(A.Grid.*B.Grid,'Size',[A.Size(1:end-order) B.Size(order+1:end)]) ;
            C.Fourier = A.Fourier || B.Fourier ;
            if isempty(A.Data) || isempty(B.Data) ; return ; end
            Ma = A.Size(1:end-order) ; 
            Na = A.Size(end-order+1:end) ; 
            Mb = B.Size(1:order) ;
            Nb = B.Size(order+1:end) ; 
            Ad = reshape(A.Data,[A.Grid.N ones(1,ndims(C.Grid)-ndims(A.Grid)) prod(Ma) prod(Na) 1 1 ]) ;
            Bd = reshape(B.Data,[B.Grid.N ones(1,ndims(C.Grid)-ndims(B.Grid)) 1 prod(Mb) prod(Nb) 1 ]) ;
            C.Data = sum(Ad.*Bd,ndims(C.Grid)+2) ; 
        end

        function x = conj(x)
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            x.Data = conj(x.Data) ;
        end

        function x = transpose(x)
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            data = x.Data ;
            switch x.Order
                case 0 % do nothing
                case 1 % column vector-> row matrix
                    x.Size = [1 x.Size] ;
                otherwise
                    perm = [2 1 3:x.Order] ;
                    data = permute(data,[1:ndims(x.Grid) ndims(x.Grid)+perm]) ;
                    x.Size = x.Size(perm) ;
            end
            x.Data = data ;
        end

        function x = ctranspose(x)
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            x = conj(transpose(x)) ;
        end

        function g = grad(x,nabla)
            if nargin<2 ; nabla = pkg.fft.fields.nabla(x.Grid) ; end
            g = tensorprod(nabla,fft(x),0) ;
        end

        function sg = symgrad(x,nabla) % faster than x = sym(grad(x)) because avoids permutations
            if nargin<2 ; nabla = pkg.fft.fields.nabla(x.Grid) ; end
            sg = pkg.fft.Field(x.Grid,'Size',[x.Size(1) x.Size(1) x.Size(2:end)],'Fourier',true) ;
            x = fft(copy(x)) ;
            if isempty(x.Data) ; return ; end
            sg.Data = .5*( nabla.Data.*reshape(x.Data,[x.Grid.N 1 x.Size]) ...
                         + reshape(nabla.Data,[nabla.Grid.N 1 nabla.Size]).*reshape(x.Data,[x.Grid.N x.Size(1) 1 x.Size(2:end)]) ) ;
        end
            
        function d = div(x,nabla)
            if nargin<2 ; nabla = pkg.fft.fields.nabla(x.Grid) ; end
            d = tensorprod(conj(nabla),fft(x),1) ;
        end

        function x = sym(x)
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            x.Data = .5*(x.Data + permute(x.Data,[1:ndims(x.Grid) ndims(x.Grid)+[2 1 3:x.Order]])) ;
            % x = 0.5.*(x+x.') ;
        end

        function t = trace(x)
            if x.Order<2 ; error('TRACE applies to tensors of order>2 !') ; end
            if x.Size(1)~=x.Size(2) ; error('TRACE applies to tensors having the two first dimensions with the same size !') ; end 
            t = pkg.fft.Field(x.Grid,'Size',x.Size(3:end),'Fourier',x.Fourier) ;
            if isa(x.Data,'gpuArray')
                t.setdata(@zeros) ;
                xd = reshape(x.Data,[prod(x.Grid.N) x.Size]) ;
                for d = 1:x.Size(1) 
                    t.Data = t.Data(:) + reshape(xd(:,d,d,:),[],1) ; 
                end 
            else
                Id = eye(x.Size(1),class(x.Data)) ;
                t.Data = sum(reshape(x.Data,prod(x.Grid.N),prod(x.Size(1:2)),prod(x.Size(3:end))).*Id(:)',2) ;
            end
        end

        function x = norm(x)
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            data = sqrt(sum(abs(reshape(x.Data,[x.Grid.N prod(x.Size) 1 1])).^2,ndims(x.Grid)+1)) ;
            x.Order = 0 ;
            x.Data = data ;
        end

        function x = inv(x)
            if x.Order~=2 ; error('INV applies to second-order tensors only !') ; end
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            x.Data = permute(pageinv(permute(x.Data,[ndims(x.Grid)+(1:2) 1:ndims(x.Grid)])),[2+(1:ndims(x.Grid)) 1 2]) ;
        end

        function x = gridfun(x,fun,dims)
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            if nargin<3 ; dims = 1:x.Grid.ndims ; end
            x.Grid = contract(x.Grid,dims) ;
            x.Data = fun(x.Data,dims) ;
        end
        function varargout = gridmean(x,varargin) ; [varargout{1:nargout}] = gridfun(x,@mean,varargin{:}) ; end
        function varargout = gridmedian(x,varargin) ; [varargout{1:nargout}] = gridfun(x,@median,varargin{:}) ; end
        function varargout = gridsum(x,varargin) ; [varargout{1:nargout}] = gridfun(x,@sum,varargin{:}) ; end

        function x = tensorfun(x,fun,dims)
            x = softcopy(x,nargout) ; % might do in-place computations if possible
            if nargin<3 ; dims = 1:x.Order ; end
            data = x.Data ; % backup
            x.Size(dims) = 1 ; % this might clear the data
            x.Data = fun(data,ndims(x.Grid)+dims) ;
        end
        function varargout = tensormean(x,varargin) ; [varargout{1:nargout}] = tensorfun(x,@mean,varargin{:}) ; end
        function varargout = tensormedian(x,varargin) ; [varargout{1:nargout}] = tensorfun(x,@median,varargin{:}) ; end
        function varargout = tensorsum(x,varargin) ; [varargout{1:nargout}] = tensorfun(x,@sum,varargin{:}) ; end
        function varargout = tensormin(x,y,varargin) 
            if nargin<2 || isempty(y)
                if x.Order==0
                    [varargout{1:nargout}] = softcopy(x,nargout) ; 
                    return ;
                else
                    y = [] ; 
                end 
            end
            [varargout{1:nargout}] = tensorfun(x,@(x,dims)min(x,y,dims),varargin{:}) ; 
        end

        function varargout = plus(varargin) ; [varargout{1:nargout}] = elemfun(@plus,varargin{:}) ; end
        function varargout = uplus(varargin) ; [varargout{1:nargout}] = elemfun(@uplus,varargin{:}) ; end
        function varargout = minus(varargin) ; [varargout{1:nargout}] = elemfun(@minus,varargin{:}) ; end
        function varargout = uminus(varargin) ; [varargout{1:nargout}] = elemfun(@uminus,varargin{:}) ; end
        function varargout = times(varargin) ; [varargout{1:nargout}] = elemfun(@times,varargin{:}) ; end
        function varargout = real(varargin) ; [varargout{1:nargout}] = elemfun(@real,varargin{:}) ; end
        function varargout = imag(varargin) ; [varargout{1:nargout}] = elemfun(@imag,varargin{:}) ; end
        function varargout = abs(varargin) ; [varargout{1:nargout}] = elemfun(@abs,varargin{:}) ; end
        function varargout = sign(varargin) ; [varargout{1:nargout}] = elemfun(@sign,varargin{:}) ; end
        function varargout = power(varargin) ; [varargout{1:nargout}] = elemfun(@power,varargin{:}) ; end

    end

end

