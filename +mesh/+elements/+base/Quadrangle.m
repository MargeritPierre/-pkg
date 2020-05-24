classdef Quadrangle < pkg.mesh.elements.AbstractElement
%QUADRANGLE 2D quad doubly-curved element with 4 nodes 

%% (MANDATORY) ELEMENT DEFINITION
% ELEMENT GEOMETRY
    properties (SetAccess = protected)
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates = [0 0 ; 1 0 ; 1 1 ; 0 1]
        % The face is itself (see constructor)
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
            N = [(1-E(:,1)).*(1-E(:,2)) E(:,1).*(1-E(:,2)) E(:,1).*E(:,2) (1-E(:,1)).*E(:,2)] ;
        end
        
        function DER = evalDerivativeAt(this,E,ORD,~)
        % Evaluate the partial derivatives of order ORD = [o1 o2 ...] (:[1 nDims]) of 
        % shape functions at given local coordinates E:[nRows nDims]
        % exemple: df(E)/(dx�dy) = DER(E,[2 1 0]) ;
            if nargin<3 ; ORD = [1 0] ; end
            if numel(ORD)~=2 ; error('Wrong derivation order argument (must be [1 nDims])') ; end
            if all(ORD==0) % No Derivative
                DER = this.evalAt(E) ;
            elseif any(ORD>1) % High-order derivatives vanish
                DER = zeros(size(E,1),this.nNodes) ;
            else % Valid derivatives
                if isequal(ORD,[1 0])
                    DER = [-(1-E(:,2)) (1-E(:,2))  E(:,2) -E(:,2)] ;
                elseif isequal(ORD,[0 1])
                    DER = [-(1-E(:,1)) -E(:,1) E(:,1) (1-E(:,1))] ;
                elseif isequal(ORD,[1 1])
                    DER = ones(size(E,1),1).*[1 -1  1 -1] ;
                end
            end
        end
    end
    
% INTEGRATION QUADRATURE
    properties
        GaussIntegrationPoints = [1/2 1/2] % [nGaussIntPts nDims]
        GaussIntegrationWeights = 1 % [nGaussIntPts 1]
    end
    
%% CONSTRUCTOR / DESTRUCTOR
    methods
        function this = Quadrangle(varargin)
        % Constructor
            this.Edges = pkg.mesh.elements.ElementTable('Types',pkg.mesh.elements.base.Bar,'Indices',[1 1 2 ; 1 2 3 ; 1 3 4 ; 1 4 1]) ;
            this.Faces = pkg.mesh.elements.ElementTable('Types',this,'Indices',[1 1 2 3 4]) ;
        end

        function delete(this)
        % Destructor
        end
    end

end

