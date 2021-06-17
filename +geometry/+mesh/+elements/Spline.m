classdef Spline < pkg.geometry.mesh.elements.AbstractElement
%Spline Bi-Spline Elements

%% (MANDATORY) ELEMENT DEFINITION
% ELEMENT GEOMETRY
    properties (SetAccess = protected)
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates
        % The list of faces [nFaces nMaxNodesByFace]
        Faces
        % The list of edges [nEdges nMaxNodesByEdge] 
        Edges
    end
    
% SHAPE FUNCTIONS AND DERIVATIVES
    methods
        function Bik = coxDeBoor(~,x,k,n)
        % Evaluate the 1D Spline functions using Cox-deBoor recursion formula
        % see https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
        % e is the local coordinates in (0,1) [nE 1]
        % k is the scalar Spline order (degree)
        % n is the scalar Spline number of knots (by default==k+1)
        % Bik is the value of the Splines at point e [nE n]
            if nargin<4 ; n = k+1 ; end 
        % Knot vector (assuming it is uniform)
            t = (0:n-1)/(n-1) ;
        % Last knot abscissa
            ii = sum(x-t>=0,2) ;
            dx = x-ti ;
        % Zero-order Splines
            Bik =  ; % [nE n+1]
        % Higher-order Splines
            for kk = 1:k
                %Bik = (n/k)*( (e-ti).*Bik(:,1:end-1) + (ti+(k+1)/n-e).*Bik(:,2:end) ) ;
                %Bik = (e-ti).*Bik(:,1:end-1) + (ti+(k+1)/(n-1)-e).*Bik(:,2:end) ; ;
                %Bik = (e-ti).*Bik ;
                %if kk<=k ; Bik = padarray(Bik,[0 1],0,'post') ; end
            end
        end
        
        % Evaluate the shape functions matrix N at local coordinates E
        % So that f(E) = N*f_nodes
        % E = [nE nDims] , N = [nE nNodes]
        function N = evalAt(this,E)
            N = ones(size(E,1),1) ;
            for d = 1:this.nDims
                Bik = this.coxDeBoor(E(:,d),this.Order(d)) ;
                N = repelem(N,1,this.Order(d)+1).*repmat(Bik,1,size(N,2)) ;
            end
        end
        
        function DER = evalDerivativeAt(this,E,ORD,~)
        % Evaluate the partial derivatives of order ORD = [o1 o2 ...] (:[1 nDims]) of 
        % shape functions at given local coordinates E:[nRows nDims]
        % exemple: df(E)/(dx²dy) = DER(E,[2 1 0]) ;
            if nargin<3 ; ORD = [1 zeros(1,this.nDims-1)] ; end
            if numel(ORD)~=this.nDims ; error('Wrong derivation order argument (must be [1 nDims])') ; end
            if all(ORD==0) % No Derivative
                DER = this.evalAt(E) ;
            elseif any(ORD>2) % High-order derivatives vanish
                DER = zeros(size(E,1),this.nNodes) ;
            else % Valid derivatives
                % See p.21 of iut.univ-lemans.fr/ydlogi/cours/mef_elas_2d.pdf
                % L = [L1(e1) L1(e2) L2(e1) L2(e2) L3(e1) L3(e2)]
                DL = [E.*(E-1)/2 1-E.^2 E.*(E+1)/2] ; % zeroth order
                DL(:,:,end+1) = [E-1/2 -2*E E+1/2] ; % first order
                if any(ORD==2) ; DL(:,:,end+1) = repmat([1 1 -2 -2 1 1],[size(E,1) 1]) ; end % second order
                % Apply the derivative
                DER = DL(:,[1 3 5 5 5 3 1 1 3],ORD(1)+1).*DL(:,[2 2 2 4 6 6 6 4 4],ORD(2)+1) ;
            end
        end
    end
    
% INTEGRATION QUADRATURE
    properties
        GaussIntegrationPoints % [nGaussIntPts nDims]
        GaussIntegrationWeights % [nGaussIntPts 1]
    end
    
%% CONSTRUCTOR / DESTRUCTOR
    properties
        Geometry = '1D' % element geometry: '1D', 'patch' or 'volume'
        Order = 3 % element order [nDims 1]
    end
    methods
        function this = Spline(geo,ord)
        % Constructor
        % Geometry
            if nargin>=1 ; this.Geometry = geo ; end
            switch this.Geometry
                case '1D' ; ndims = 1 ;
                case 'patch' ; ndims = 2 ;
                case 'volume' ; ndims = 3 ;
            end
        % Order 
            if nargin>=2 ; this.Order = ord ; end
            this.Order = this.Order(:).*ones(ndims,1) ;
        % Node coordinates
            X = arrayfun(@linspace,zeros(ndims,1),ones(ndims,1),this.Order+1,'UniformOutput',false) ;
            [X{:}] = ndgrid(X{:}) ;
            this.NodeLocalCoordinates = reshape(cat(ndims+1,X{:}),[],ndims) ;
        end

        function delete(this)
        % Destructor
        end
    end

end

%% TEST FUNCTIONS
function test

%% CREATE A 1D BI-SPLINE ELEMENT AND PLOT ITS SHAPE FUNCTIONS
elmt = pkg.geometry.mesh.elements.Spline('1D',2) ;
E = linspace(0,1,100)' ;
N = elmt.evalAt(E) ;
clf ; plot(E,N) ;



end
