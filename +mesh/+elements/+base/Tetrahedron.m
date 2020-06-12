classdef Tetrahedron < pkg.mesh.elements.base.BaseElement
%TETRAHEDRON 3D simplex with 4 nodes

%% (MANDATORY) ELEMENT DEFINITION
% ELEMENT GEOMETRY
    properties (SetAccess = protected)
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates = [0 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1]
        % The list of faces [nFaces 3]
        Faces
        % The list of edges [nEdges 2] 
        Edges
    end
    
% SHAPE FUNCTIONS AND DERIVATIVES
    methods
        % Evaluate the shape functions matrix N at local coordinates E
        % So that f(E) = N*f_nodes
        % E = [nE nDims] , N = [nE nNodes]
        function N = evalAt(~,E)
            N = [1-E(:,1)-E(:,2)-E(:,3) E(:,1) E(:,2) E(:,3)] ;
        end
        
        function DER = evalDerivativeAt(this,E,ORD,~)
        % Evaluate the partial derivatives of order ORD = [o1 o2 ...] (:[1 nDims]) of 
        % shape functions at given local coordinates E:[nRows nDims]
        % exemple: df(E)/(dx²dy) = DER(E,[2 1 0]) ;
            if nargin<3 ; ORD = [1 0 0] ; end
            if numel(ORD)~=3 ; error('Wrong derivation order argument (must be [1 nDims])') ; end
            if all(ORD==0) % No Derivative
                DER = this.evalAt(E) ;
            elseif sum(ORD)>1 % High-order derivatives vanish
                DER = zeros(size(E,1),this.nNodes) ;
            else % Valid derivatives
                if isequal(ORD,[1 0 0])
                    DER = ones(size(E,1),1).*[-1 1 0 0] ;
                elseif isequal(ORD,[0 1 0])
                    DER = ones(size(E,1),1).*[-1 0 1 0] ;
                elseif isequal(ORD,[0 0 1])
                    DER = ones(size(E,1),1).*[-1 0 0 1] ;
                end
            end
        end
    end
    
% INTEGRATION QUADRATURE
    properties
        GaussIntegrationPoints = [1/3 1/3 1/3] % [nGaussIntPts nDims]
        GaussIntegrationWeights = 1/6 % [nGaussIntPts 1]
    end
    
%% CONSTRUCTOR / DESTRUCTOR
    methods
        function this = Tetrahedron(varargin)
        % Constructor
            this.Edges = pkg.mesh.elements.ElementTable('Types',pkg.mesh.elements.base.Bar...
                                                            ,'Indices',[1 1 2 ; 1 2 3 ; 1 3 1 ; 1 1 4 ; 1 2 4 ; 1 3 4]) ;
            this.Faces = pkg.mesh.elements.ElementTable('Types',pkg.mesh.elements.base.Triangle...
                                                            ,'Indices',[1 1 3 2 ; 1 1 2 4 ; 1 2 3 4 ; 1 3 1 4]) ;
        end

        function delete(this)
        % Destructor
        end
    end

end

