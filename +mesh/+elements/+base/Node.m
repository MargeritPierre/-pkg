classdef Node < pkg.mesh.elements.AbstractElement
%NODE 0D element with one node

%% (MANDATORY) ELEMENT DEFINITION
% ELEMENT GEOMETRY
    properties (SetAccess = protected)
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates = [0]
        % No faces
        Faces = pkg.mesh.elements.ElementTable
        % No Edges
        Edges = pkg.mesh.elements.ElementTable 
    end
    
% SHAPE FUNCTIONS AND DERIVATIVES
    methods
        % Evaluate the shape functions matrix N at local coordinates E
        % So that f(E) = N*f_nodes
        % E = [nE nDims==0] , N = [nE nNodes==1]
        function N = evalAt(~,E)
            N = ones(size(E,1),1) ;
        end
        
        function DER = evalDerivativeAt(this,E,ORD,~)
        % Evaluate the partial derivatives of order ORD = [o1 o2 ...] (:[1 nDims]) of 
        % shape functions at given local coordinates E:[nRows nDims]
        % exemple: df(E)/(dx²dy) = DER(E,[2 1 0]) ;
            if nargin<3 ; ORD = 1; end
            if numel(ORD)~=1 ; error('Wrong derivation order argument (must be scalar)') ; end
            if ORD==0 % No Derivative
                DER = this.evalAt(E) ;
            else % High-order derivatives vanish
                DER = zeros(size(E,1),2) ;
            end
        end
    end
    
% INTEGRATION QUADRATURE
    properties
        GaussIntegrationPoints = 0 % [nGaussIntPts nDims]
        GaussIntegrationWeights = 0 % [nGaussIntPts 1]
    end
    
%% CONSTRUCTOR / DESTRUCTOR
    methods
        function this = Node(varargin)
        % Constructor
        end

        function delete(this)
        % Destructor
        end
    end

end

