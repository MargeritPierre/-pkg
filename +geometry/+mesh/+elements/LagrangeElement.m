classdef LagrangeElement < pkg.geometry.mesh.elements.AbstractElement
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
        Type
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates
        % The list of faces [nFaces nMaxNodesByFace]
        Faces = pkg.geometry.mesh.elements.ElementTable ;
        % The list of edges [nEdges 2] 
        Edges = pkg.geometry.mesh.elements.ElementTable ;
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
                    case 'tri' % 2D simplex triangle
                        this.createSimplexElement(2,order) ;
                    case 'quad' % 2D quad
                        this.createTensorProdElement(2,order) ;
                    case 'tet' % 3D tetrahedron
                        this.createSimplexElement(3,order) ;
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
        
        function edges = faceEdges(~,faces,edges)
        % Return the node indices of unique edges related to faces
            nFaces = size(faces,1) ; nEdges = size(edges,1) ; nNd = size(edges,2) ;
            edges = sub2ind(size(faces),repmat((1:nFaces),[numel(edges) 1]),repmat(edges(:),[1 nFaces])) ;
            edges = reshape(faces(edges),[nEdges nNd nFaces]) ;
            edges = reshape(permute(edges,[1 3 2]),[nFaces*nEdges nNd]) ;
            [~,ue] = unique(sort(edges,2),'rows') ;
            edges = edges(ue,:) ;
        end
        
        function createTensorProdElement(this,NDIMS,ORDER)
        % Set the element properties corresponding to a tensor product element 
        % Used for geometries 1D, quad, hex, ...
            this.Type = 'tensorprod' ;
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
        % ELEMENT EDGES AND FACES
            nNd = max(2,this.Order+1) ; % number of nodes by dimension
            IND = reshape(1:nNd^this.nDims,[nNd*ones(1,this.nDims) 1]) ; % Array of node indices
            % 1D elements
                if this.nDims==1
                    edgeElemType = this ; % The element is 1D, so its only edge is itself
                    edgeNodIdx = IND(:)' ; 
                else
                    edgeElemType = pkg.geometry.mesh.elements.LagrangeElement('1D',this.Order) ; % Edges: 1D element of the same order
                    e1 = IND(:,1,1) ; e2 = IND(end,:,1) ; e3 = flip(IND(:,end,1)) ; e4 = flip(IND(1,:,1)) ;
                    edgeNodIdx = [e1(:) e2(:) e3(:) e4(:)]' ;
                end
            % 2D elements
                if this.nDims==2
                    faceElemType = this ; % The element is 2D, so its only face is itself 
                    faceNodIdx = 1:nNd^2 ;
            % 3D elements 
                elseif this.nDims>2
                    faceElemType = pkg.geometry.mesh.elements.LagrangeElement('quad',this.Order) ; % Faces: 2D quad of the same order
                    % Faces
                        f1 = permute(IND(1,:,:),[3 2 1]) ; 
                        f2 = IND(end,:,:) ; 
                        f3 = IND(:,1,:) ; 
                        f4 = permute(IND(:,end,:),[3 1 2]) ; 
                        f5 = IND(:,:,1)' ; 
                        f6 = IND(:,:,end) ;
                        faceNodIdx = [f1(:) f2(:) f3(:) f4(:) f5(:) f6(:)]' ;
                    % Edges
                        edgeNodIdx = faceEdges(this,faceNodIdx,edgeNodIdx) ;
                end
            % Assign
                if this.nDims>1 ; this.Faces = pkg.geometry.mesh.elements.ElementTable('Types',faceElemType,'NodeIdx',faceNodIdx) ; end
                this.Edges = pkg.geometry.mesh.elements.ElementTable('Types',edgeElemType,'NodeIdx',edgeNodIdx) ;
        % GAUSS QUADRATURE POINTS
            this.setGaussIntegration() ;
        end
        
        function createSimplexElement(this,NDIMS,ORDER)
        % Set the element properties corresponding to a simplex element 
        % Used for geometries tri, tet, ...
            this.Type = 'simplex' ;
        % CONSTRUCT THE BASE POLYNOMS
            if ORDER==0 % ZERO-order elements
                this.PolyExp = zeros(1,NDIMS) ;
                this.NodeLocalCoordinates = [zeros(1,NDIMS) ; eye(NDIMS)] ;
            else % Higher-order elements
                % Polynomial coefficients
                    P = arrayfun(@(dim)0:ORDER,1:NDIMS,'UniformOutput',false) ;
                    if NDIMS>1 ; [P{:}] = ndgrid(P{:}) ; end
                    PP = cat(NDIMS+1,P{:}) ;
                    PP = reshape(PP,(ORDER+1)^NDIMS,NDIMS) ;
                    this.PolyExp = PP(sum(PP,2)<=ORDER,:) ;
                % Node coordinates
                    this.NodeLocalCoordinates = this.PolyExp/ORDER ;
            end
            this.buildPolynomialBasis() ;
        % ELEMENT EDGES AND FACES
            nNd = max(2,this.Order+1) ; % number of nodes by dimension
            % Edges
                edgeElemType = pkg.geometry.mesh.elements.LagrangeElement('1D',this.Order) ; % Edges: 1D element of the same order
                e1 = 1:nNd ;
                e2 = cumsum(nNd:-1:1) ;
                e3 = flip([1 e2(1:end-2)+1 e2(end)]) ;
                edgeNodIdx = [e1 ; e2 ; e3] ;
            % 2D Simplices
                if this.nDims==2
                    faceElemType = this ; % The element is 2D, so its only face is itself 
                    faceNodIdx = 1:this.nNodes ;
            % 3D elements 
                elseif this.nDims>2
                    IND = reshape(PP,[nNd*ones(1,this.nDims) this.nDims]) ;
                    IND = double(sum(IND,this.nDims+1)<nNd) ;
                    IND(IND~=0) = 1:this.nNodes ;
                    first = [1 zeros(1,nNd-1)] ; 
                    p2IND = permute(IND,[1 3 2]) ;
                    p3IND = permute(IND,[2 1 3]) ;
                    % Faces
                        faceElemType = pkg.geometry.mesh.elements.LagrangeElement('tri',this.Order) ; % Faces: 2D triangles of the same order
                        f1 = IND((IND.*first)~=0) ; %find(this.NodeLocalCoordinates(:,1)==0) ;
                        f2 = p2IND((p2IND.*first')~=0) ; %find(this.NodeLocalCoordinates(:,2)==0) ;
                        f3 = p3IND((p3IND.*reshape(first,1,1,[]))~=0) ; %find(this.NodeLocalCoordinates(:,3)==0) ;
                        f4 = find(abs(sum(this.NodeLocalCoordinates,2)-1)<eps) ;
                        faceNodIdx = [f1(:) f2(:) f3(:) f4(:)]' ;
                    % Edges
                        edgeNodIdx = faceEdges(this,faceNodIdx,edgeNodIdx) ;
                end
            % Assign
                this.Faces = pkg.geometry.mesh.elements.ElementTable('Types',faceElemType,'NodeIdx',faceNodIdx) ;
                this.Edges = pkg.geometry.mesh.elements.ElementTable('Types',edgeElemType,'NodeIdx',edgeNodIdx) ;
        % GAUSS QUADRATURE POINTS
            this.setGaussIntegration() ;
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
        % be given to evaluate polynom derivatives
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
    
    
%% INTEGRATION QUADRATURE
    properties
        GaussIntegrationPoints % [nGaussIntPts nDims]
        GaussIntegrationWeights % [nGaussIntPts 1]
    end
    methods
        function setGaussIntegration(this,order)
        % Compute the integration scheme
        % see www.code-aster.org › doc › man_r › r3.01.01.pdf
            if nargin<2 ; order = this.Order+1 ; end
            switch this.Type
                case 'tensorprod'
                % 1D Case
                    switch order
                        case 0 ; this.setGaussIntegration(1) ; return ; 
                        case 1 ; GP = 1/2 ; W = 1 ;
                        case 2
                            a = 1/sqrt(3) ;
                            GP = ([-a;a]+1)/2 ;
                            W = [1;1]/2 ;
                        case 3
                            a = sqrt(3/5) ;
                            GP = ([-a;0;a]+1)/2 ;
                            W = [5;8;5]/18 ;
                        otherwise
                            warning('Element integration scheme limited to order 3.') ;
                            this.setGaussIntegration(3) ;
                            return ; 
                    end
                % Multidim grid
                    for dd = 2:this.nDims
                        GP = [repmat(GP,order,1) repelem(GP(:,end),order,1)] ;
                        W = repmat(W(:),order,1).*repelem(W(:),order,1) ;
                    end
                case 'simplex'
                    switch order
                        case 0 ; this.setGaussIntegration(1) ; return ; 
                        case 1
                            GP = this.centroid ;
                            W = (1/factorial(this.nDims)) ;
                        case 2
                            switch this.nDims
                                case 2 % FPG3
                                    GP = [1 1 ; 4 1 ; 1 4]/6 ;
                                    W = [1;1;1]/6 ;
                                case 3 % FPG4
                                    a = (5-sqrt(5))/30 ;
                                    b = (5+3*sqrt(5))/20 ;
                                    GP = [a a a ; a a b ; a b a ; b a a] ;
                                    W = [1;1;1;1]/24 ;
                            end
                        case 3
                            switch this.nDims
                                case 2 % FPG6
                                    a  = 0.445948490915965 ; 
                                    b  = 0.091576213509771 ;
                                    GP = [b b ; 1-2*b b ; b 1-2*b ; a 1-2*a ; a a ; 1-2*a a] ;
                                    P1 = 0.11169079483905 ;
                                    P2 = 0.0549758718227661 ; 
                                    W = [P2;P2;P2;P1;P1;P1] ;
                                case 3 % FPG5
                                    a = 0.25 ;
                                    b = 1/6 ;
                                    c = 0.5 ;
                                    GP = [a a a ; b b b ; b b c ; b c b ; c b b] ;
                                    W = [-2/15;[1;1;1;1]*3/40] ;
                            end
                        otherwise
                            warning('Element integration scheme limited to order 3.') ;
                            this.setGaussIntegration(3) ;
                            return ; 
                    end
            end
        % Set
            this.GaussIntegrationPoints = GP ;
            this.GaussIntegrationWeights = W ;
        end
    end

end

