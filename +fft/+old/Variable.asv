classdef Variable < pkg.fft.Space
% A tensor field variable (Space+Data) defined on a regular grid

    properties 
        Data
    end
    
    methods
        function this = Variable(varargin)
        % Constructor
        % this = pkg.fft.Variable(space) ; % transparent
        % this = pkg.fft.Variable(data) ; % conversion of some data to a "constant" Variable
        % this = pkg.fft.Variable(grid) ;
        % this = pkg.fft.Variable(...,Name,Value) ;
            this = this@pkg.fft.Space(varargin{:}) ;
            if nargin==0 ; return ; end
            if isnumeric(varargin{1}) ; this.Data = varargin{1} ; end
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

        function pl = plot(varargin)
            if nargin>1 && isa(varargin{2},'pkg.fft.Variable') % plot(x,variable,...)
                x = varargin{1} ;
                var = varargin{2} ;
                varargin = varargin(3:end) ;
            else % plot(variable,...)
                var = varargin{1} ;
                x = var.Grid.coordinates ;
                varargin = varargin(2:end) ;
            end
            if var.Fourier
                var = copy(var) ; % to not modify the plotted variable data
                x = var.Grid.fouriercoordinates ; % change the plotted coordinates to Fourier
                for d = 1:var.Grid.ndims
                    var.Data = fftshift(var.Data,d) ;
                    x.Data = fftshift(x.Data,d) ;
                end
            end
            pl = plot(var.Grid,varargin{:}) ;
            pl.FaceColor = 'interp' ;
            pl.EdgeColor = 'none' ;
            pl.XData = gather(x.Data(:,:,1)) ;
            pl.YData = gather(x.Data(:,:,2)) ;
            if var.Order==0
                pl.CData = gather(var.Data) ;
            else
                pl.CData = gather(tensornorm(var).Data) ;
            end
        end
    
    end


    methods

        function x = datacast(x,fcn)
            x = datacast@pkg.fft.Space(x,fcn) ;
            x.Data = fcn(x.Data) ;
        end

        function x = fouriertransform(x,inverse)
            x = fouriertransform@pkg.fft.Space(x,inverse) ;
            if ~isempty(x.Data)
                for d = 1:x.Grid.ndims
                    if inverse
                        x.Data = ifft(x.Data,x.Grid.N(d),d) ;
                    else
                        x.Data = fft(x.Data,x.Grid.N(d),d) ;
                    end
                end
            end
        end

        function C = arrayfun(fun,varargin)
            C = pkg.fft.Variable(arrayfun@pkg.fft.Space(fun,varargin{:})) ;
        % convert all inputs to spaces
            isvar = cellfun(@(arg)isa(arg,'pkg.fft.Variable'),varargin) ;
            if ~all(isvar)
                varargin = cellfun(@pkg.fft.Variable,varargin) ;
            else
                varargin = [varargin{:}] ;
            end
        % Apply function on data
            if numel(unique([varargin.Grid]))==1 % all fields are on the same grid
                data = fun(varargin.Data) ; 
            else
                data = cellfun(@(data,grid,sz)reshape(data,[grid.N ones(1,ndims(C.Grid)-ndims(grid)) sz]),{varargin.Data},{varargin.Grid},{varargin.Size},'uni',false) ;
                data = fun(data{:}) ;
            end
            C.Data = data ;
        end

        function C = tensorprod(A,B,order)
        % TENSORPROD of two gridded tensor field variables at a given order
            C = pkg.fft.Variable(tensorprod@pkg.fft.Space(A,B,order)) ;
            if isempty(A.Data) || isempty(B.Data) ; return ; end
            Ma = A.Size(1:end-order) ; 
            Na = A.Size(end-order+1:end) ; 
            Mb = B.Size(1:order) ;
            Nb = B.Size(order+1:end) ; 
            Ad = reshape(A.Data,[A.Grid.N ones(1,ndims(C.Grid)-ndims(A.Grid)) prod(Ma) prod(Na) 1 1 ]) ;
            Bd = reshape(B.Data,[B.Grid.N ones(1,ndims(C.Grid)-ndims(B.Grid)) 1 prod(Mb) prod(Nb) 1 ]) ;
            C.Data = sum(Ad.*Bd,ndims(C.Grid)+2) ; 
        end

        function x = transpose(x)
            data = x.Data ;
            x = pkg.fft.Variable(transpose@pkg.fft.Space(x)) ;
            switch x.Order
                case 0 % do nothing
                case 1 % column vector-> row matrix
                otherwise
                    perm = [2 1 3:x.Order] ;
                    data = permute(data,[1:ndims(x.Grid) ndims(x.Grid)+perm]) ;
            end
            x.Data = data ;
        end

        function t = trace(x)
            t = pkg.fft.Variable(trace@pkg.fft.Space(x)) ;
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

        function y = gridfun(x,fun,dims)
            y = pkg.fft.Variable(gridfun@pkg.fft.Space(x)) ;
            y.Data = fun(x.Data,dims) ;
        end

        function y = tensorfun(x,fun,dims)
            y = pkg.fft.Variable(tensorfun@pkg.fft.Space(x)) ;
            y.Data = fun(data,ndims(x.Grid)+dims) ;
        end

        function y = tensornorm(x)
            y = pkg.fft.Variable(tensornorm@pkg.fft.Space(x)) ;
            y.Data = sqrt(sum(abs(reshape(x.Data,[x.Grid.N prod(x.Size) 1 1])).^2,ndims(x.Grid)+1)) ;
        end

        function y = tensorinv(x)
            y = pkg.fft.Variable(tensorinv@pkg.fft.Space(x)) ;
            y.Data = permute(pageinv(permute(x.Data,[ndims(x.Grid)+(1:2) 1:ndims(x.Grid)])),[2+(1:ndims(x.Grid)) 1 2]) ;
        end

    end

end

