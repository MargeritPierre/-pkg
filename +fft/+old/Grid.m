classdef Grid < matlab.mixin.Copyable
% GRID A Uniform grid
    
    properties
        N(1,:) {mustBeInteger, mustBeNonnegative} ; % number of elements
        x0(:,1) % grid origin
        dx(:,:) % grid steps or unit vectors
    end
    
    methods
        function this = Grid(N,varargin)
        % Constructor
            this.N = N ;
            for arg = 1:2:numel(varargin)
                this.(varargin{arg}) = varargin{arg+1} ;
            end
            if isempty(this.x0)
                D = max(numel(this.N),size(this.dx,2)) ;
                this.x0 = zeros([D 1],class(this.dx)) ;
            end
            this.N = this.N + zeros([1 size(this.x0,1)],class(this.N)) ;
            if isempty(this.dx)
                this.dx = eye(numel(this.N),class(this.x0)) ;
            end
            if isvector(this.dx)
                this.dx = diag(this.dx) ;
            end
        end

        function D = ndims(this) ; D = uint32(numel(this.N)) ; end

        function idx = data_indices(this)
            % idx = arrayfun(@colon,0*this.N,this.N-1,'uni',false) ; % "uni,false" does not support gpuArrays
            idx = cellfun(@colon,num2cell(0*this.N,1),num2cell(this.N-1,1),'uni',false) ;
            [idx{:}] = ndgrid(idx{:}) ;
            idx = cellfun(@(i)i(:)',idx,'uni',false) ;
            idx = cat(1,idx{:}) ;
            idx = reshape(idx,[this.ndims this.N]) ;
        end

        function x = coordinates(this)
            x = pkg.fft.old.Field(this,'Order',1) ;
            idx = this.data_indices() ;
            x.Data = this.x0 + reshape(this.dx*idx(:,:),size(idx)) ;
        end

        function this = cast(this,castfcn)
            this = copy(this) ;
            this.N = castfcn(this.N) ;
            this.dx = castfcn(this.dx) ;
            this.x0 = castfcn(this.x0) ;
        end
        function this = double(this) ; this = cast(this,@double) ; end
        function this = single(this) ; this = cast(this,@single) ; end
        function this = uint32(this) ; this = cast(this,@uint32) ; end
        function this = gpuArray(this) ; this = cast(this,@gpuArray) ; end

        function this = fourier(this)
            this = copy(this) ;
            this.dx = 2*pi*((this.dx.*this.N)'\eye(this.ndims)) ;
            this.x0 = -this.dx*floor(this.N(:)/2) ;
        end

        function xi = fouriercoordinates(this)
            fouriergrid = this.fourier() ;
            xi = fouriergrid.coordinates() ;
            for d = 1:this.ndims
                xi.Data = ifftshift(xi.Data,d+xi.Order) ;
            end
        end

        function grid = contract(this,dims)
            grid = copy(this) ;
            grid.N(dims) = 1 ;
            grid.dx = grid.dx.*(this.N./grid.N) ;
        end

        function Gc = times(Ga,Gb)
            if Ga==Gb ; Gc = Ga ; return ; end
            if Ga.ndims==0 || all(Ga.N==1) ; Gc = Gb ; return ; end
            if Gb.ndims==0 || all(Gb.N==1) ; Gc = Ga ; return ; end
            Na = [Ga.N zeros(1,max(0,Gb.ndims-Ga.ndims))] ;
            Nb = [Gb.N zeros(1,max(0,Ga.ndims-Gb.ndims))] ;
            valid = Na==Nb | Na<=1 | Nb<=1 ;
            if ~all(valid) ; error('Incompatible grid dimensions for tensor product !') ; end
            Gc = pkg.fft.old.Grid(max(Na,Nb)) ;
            dxA = padarray(Ga.dx,max(0,size(Gb.dx)-size(Ga.dx)),0,'post') ;
            dxB = padarray(Gb.dx,max(0,size(Ga.dx)-size(Gb.dx)),0,'post') ;
            Gc.dx = dxA.*(Na>1) + dxB.*(Nb>1) ;
            Gc.dx(:,Na==Nb) = .5*(dxA(:,Na==Nb) + dxB(:,Na==Nb)) ;
            x0a = padarray(Ga.x0,max(0,size(Gb.x0)-size(Ga.x0)),0,'post') ;
            x0b = padarray(Gb.x0,max(0,size(Ga.x0)-size(Gb.x0)),0,'post') ;
            Gc.x0 = x0a + x0b ;
            % if Ga.Fourier~=Gb.Fourier ; warning('Mixing two grids that are not on different spaces !') ; end
            % Gc.Fourier = Ga.Fourier || Gb.Fourier ;
        end

        function pl = plot(this,varargin)
            x = this.coordinates.Data ;
            % if this.Fourier
            %     for d = 1:this.ndims
            %         x = fftshift(x,d+1) ;
            %     end
            % end
            switch this.ndims
                case 1
                    pl = plot(x(:),0*x(:),'.-',varargin{:}) ;
                case 2
                    pl = surf(squeeze(x(1,:,:)),squeeze(x(2,:,:)),zeros(this.N),NaN(this.N),varargin{:}) ;
                case 3
                otherwise
            end
        end
    end

end

