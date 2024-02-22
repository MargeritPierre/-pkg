classdef GridMesh < pkg.geometry.mesh.Mesh
%GRIDMESH Construct a grid mesh

%% GRID MESH CONSTRUCTOR
properties
    Ne(1,:) double = [] ; % number of elements
    Domain(2,:) double = [] ; % grid domain
end
methods
    function this = GridMesh(varargin)
    % Constructor
    % mesh = GridMesh(Ne) with Ne [1 nCoord] INT (number of elements)
    % mesh = GridMesh(L) with L [1 nCoord] DOUBLE
    % mesh = GridMesh(domain) with domain [2 nCoord] DOUBLE ([xmin ; xmax]
    % mesh = GridMesh(...,Ne) with Ne scalar or [1 nCoord] INT
    % mesh = GridMesh(...,dx) with dx scalar or [1 nCoord] DOUBLE
    % mesh = GridMesh(...,varargin) to set other mesh properties
        defN = 1000 ; % default total number of elements
    % Process first argument
        if nargin>0
            if isinteger(varargin{1}) % mesh = GridMesh(N)
                Ne = double(varargin{1}) ;
                domain = Ne.*[0;1] ;
            else
                if size(varargin{1},1)==1 % mesh = GridMesh(L,..)
                    L = varargin{1} ;
                    domain = L.*[0;1] ;
                else % mesh = GridMesh(domain,..)
                    domain = varargin{1} ;
                    L = range(domain,1) ;
                end
                Ne = ceil(L.*(defN/prod(L)).^(1/numel(L))) ;
            end
            varargin(1) = [] ;
        end
    % Second argument ?
        if nargin>1
            if isinteger(varargin{1}) % mesh = GridMesh(...,N)
                Ne = double(varargin{1}) ;
            else % mesh = GridMesh(...,dx)
                dx = varargin{1} ;
                Ne = round(range(domain,1)./dx) ;
            end
            varargin(1) = [] ;
        end
    % List of nodes and elements
        if nargin==0
            X = [] ; Elems = [] ;
        else
        % Format input
            domain = domain + 0*Ne ;
            Ne = 0*domain(1,:) + Ne ;
            nCoord = numel(Ne) ; % number of grid coordinates
            n = Ne+1 ; % number of nodes by dimension
        % Grid coordinates
            X = arrayfun(@linspace,domain(1,:),domain(2,:),n,'UniformOutput',false) ;
            [X{:}] = ndgrid(X{:}) ;
            X = cat(nCoord+1,X{:}) ;
            X = reshape(X,[],nCoord) ;
        % Element indices
            Elems = 1 ;
            for cc = find(Ne>0)
                ic = (0:Ne(cc)-1)' + [0 1] ;
                Elems = repmat(Elems,Ne(cc),2) + prod(n(1:cc-1))*repelem(ic,size(Elems,1),size(Elems,2)) ;
                if sum(Ne(1:cc)>0)==2
                    Elems = Elems(:,[1 2 4 3]) ;
                end
            end
        end
    % Mesh
        this = this@pkg.geometry.mesh.Mesh('Nodes',X,'Elems',Elems,varargin{:}) ;
        this.Ne = Ne ;
        this.Domain = domain ;
    end
end
    
end

