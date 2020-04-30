classdef LagrangeElement < pkg.mesh.elements.AbstractElement
%LAGRANGEELEMENT Implements the general Lagrange elements
% 
% elmt = LAGRANGEELEMENT(geometry,order)
% Creates a Lagrange element wit the given inputs:
% geometry: class of the element [string] (mandatory):
%   - '1D': 1D element
%   - 'tri': 2D simplex triangle
%   - 'quad': 2D quad
%   - 'tet': 3D tetrahedron
%   - 'hex': 3D hexahedron
%   - ...
% order: interpolation order [integer] (optionnal, default = 1)

%% ELEMENT GEOMETRY
    properties (SetAccess = protected)
        % Creation arguments
        Geometry
        Order
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates
        % The list of faces [nFaces nMaxNodesByFace]
        Faces ;
        % The list of edges [nEdges 2] 
        Edges ;
    end
    
%% POLYNOMIAL SHAPE FUNCTIONS
    properties (SetAccess = protected)
        % http://femwiki.wikidot.com/elements:lagrange-elements
        % Ni(x) = A_{i,n,d}*x_d^n = a_{i,m}*P_m(x)
        % i = 1,..,nNodes
        % n = 0,..,order
        % d = 1,..,nDims
        % m = 1,...,nPolys
        % Except for order==0, nPolys=nNodes
        % Base Polynomial EXPONENTS [nPolys nDims]
        PolyExp
        % Polynomial COEFFICIENTS of shape functions [nPolys nNodes]
        PolyCoeffs
    end
    methods
        % Number of polynoms in the basis
        function N = nPolys(this) ; N = size(this.PolyExp,1) ; end
    end
    
%% INTEGRATION QUADRATURE
    properties
        GaussIntegrationPoints % [nGaussIntPts nDims]
        GaussIntegrationWeights % [nGaussIntPts 1]
    end
    
%% CONSTRUCTOR / DESTRUCTOR
    methods
        function this = LagrangeElement(geometry,order)
        % Constructor of the class
            if nargin<2 ; order = 1 ; end
            % Record Inputs
                this.Geometry = geometry ;
                this.Order = order ;
            % Geometry-dependent properties
                switch geometry
                    case '1D' % 1D element
                        this.createTensorProdElement(1,order) ;
                        this.Edges = (1:this.nNodes-1)' + [0 1] ;
                    case 'tri' % 2D simplex triangle
                    case 'quad' % 2D quad
                        this.createTensorProdElement(2,order) ;
                    case 'tet' % 3D tetrahedron
                    case 'hex' % 3D hexahedron
                        this.createTensorProdElement(3,order) ;
                    otherwise
                        error('Unknown geometry argument.') ;
                end
        end
        
        function buildPolynomialBasis(this)
        % Compute the shape function coefficients
            if this.Order==0 % Zero-order elements: mean over the nodes
                this.PolyCoeffs = ones(1,this.nNodes)/this.nNodes ;
            else % General case 
                Pe = this.evalPolynomsAt(this.NodeLocalCoordinates) ; % [nNodes nPolys]
                this.PolyCoeffs = Pe\eye(this.nNodes) ; % [nPolys nNodes]
            end
        end
        
        function createTensorProdElement(this,NDIMS,ORDER)
        % Set the element properties corresponding to a tensor product element 
        % Used for geometries 1D, quad, hex, ...
            % CONSTRUCT THE BASE POLYNOMS
                if ORDER==0 % ZERO-order elements
                    this.PolyExp = zeros(1,NDIMS) ;
                    this.NodeLocalCoordinates = [zeros(1,NDIMS) ; eye(NDIMS) ; ones(1,NDIMS)] ;
                else
                    % Polynom Exponents
                        P = arrayfun(@(dim)0:ORDER,1:NDIMS,'UniformOutput',false) ;
                        if NDIMS>1 ; [P{:}] = ndgrid(P{:}) ; end
                        PP = cat(NDIMS+1,P{:}) ;
                        this.PolyExp = reshape(PP,(ORDER+1)^NDIMS,NDIMS) ;
                    % Node Coordinates
                        this.NodeLocalCoordinates = this.PolyExp/ORDER ;
                end
                this.buildPolynomialBasis() ;
            % ELEMENT CONNECTIVITIES
                nNd = max(2,this.Order+1) ; % number of nodes by dimension
                % Uses subsref to find edges, faces, etc
                    IND = reshape(1:this.nNodes,nNd*ones(1,this.nDims)) ;
                    S = [] ; S.type = '()' ; S.subs = cell(1,this.nDims) ;
                % Edges
                    this.Edges = [] ;
                    for dim = 1:this.nDims
                        for ii = [1 nNd]
                            [S.subs{:}] = deal(ii) ;
                            S.subs{dim} = ':' ;
                            iii = reshape(subsref(IND,S),nNd,1) ;
                            this.Edges = [this.Edges ; iii(1:end-1) iii(2:end)] ;
                        end
                    end
                % Faces
                    this.Faces = [] ;
%                     if this.nDims==2
%                         for dim = 1:this.nDims
%                             step = prod(nNd*(1:dim-1)) ;
%                             this.Edges = [this.Edges ; step*((0:nNd-2)' + [0 1]) + 1 ] ; % First line
%                             this.Edges = [this.Edges ; (this.nNodes:-step:this.nNodes-(nNd-2)*step)' - step*(0:1) ] ; % Last line
%                         end
%                     elseif this.nDims>2
%                     end
        end
    end
    
%% SHAPE FUNCTIONS AND DERIVATIVES
    methods
        % Evaluate the shape functions matrix N at local coordinates E
        % So that f(E) = N*f_nodes
        % E = [nE nDims] , N = [nE nNodes]
        function N = evalAt(this,E)
            P = this.evalPolynomsAt(E) ;
            N = P*this.PolyCoeffs ;
        end
        
        function DER = evalDerivativeAt(this,E,ORD,~)
        % Evaluate the partial derivatives of order ORD = [o1 o2 ...] (:[1 nDims]) of 
        % shape functions at given local coordinates E:[nRows nDims]
        % exemple: df(E)/(dx²dy) = DER(E,[2 1 0]) ;
            if nargin<3 ; ORD = [1 zeros(1,this.nDims-1)] ; end
            if numel(ORD)~=this.nDims ; error('Wrong derivation order argument (must be [1 nDims])') ; end
            if all(ORD==0) % No Derivative
                DER = this.evalAt(E) ;
            elseif any(ORD>this.Order) % High-order derivatives vanish
                DER = zeros(size(E,1),this.nNodes) ;
            else % Valid derivatives
                P = this.evalPolynomsAt(E,ORD) ;
                DER = P*this.PolyCoeffs ;
            end
        end
        
        function P = evalPolynomsAt(this,E,DERIV)
        % Evaluate base polynoms P at given local coordinates E
        % An optional derivation argument DERIV = [d1 d2 ...] [1 nDims] can
        % be gicven to evaluate polynom derivatives
        % E: [nRows nDims]
        % P: [nRows nPolys]
            if nargin<3 ; DERIV = zeros(1,this.nDims) ; end
            E = reshape(E,[size(E,1) 1 this.nDims]) ; % [nRows 1 nDims]
            N = reshape(this.PolyExp,[1 this.nPolys this.nDims]) ; % [1 nPolys nDims]
            if all(DERIV==0) % No derivatives, P_{i,j} = x{i,d}.^n{j,d}
                P = prod(E.^N,3) ; % [nRows nPolys]
            else % With derivatives
                DERIV = reshape(DERIV,[1 1 this.nDims]) ;
                % New exponents
                    Nm = N-DERIV ;
                    valid = Nm>=0 ; 
                    Nm(~valid) = 0 ;
                % Coefficients
                    A = zeros(size(Nm)) ;
                    A(valid) = factorial(N(valid))./factorial(Nm(valid)) ;
                % Compute
                    P = prod(A.*(E.^Nm),3) ; % [nRows nPolys]
            end
        end
    end

end

