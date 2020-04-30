classdef EmptyElement < pkg.mesh.elements.AbstractElement
%EMPTYELEMENT Empty Element

%% (MANDATORY) ELEMENT DEFINITION
% ELEMENT GEOMETRY
    properties (SetAccess = protected)
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates = []
        % The list of faces [nFaces nMaxNodesByFace]
        Faces = [] % nineth node is included
        % The list of edges [nEdges 2] 
        Edges = [] 
    end
    
% SHAPE FUNCTIONS AND DERIVATIVES
    methods
        % Evaluate the shape functions matrix N at local coordinates E
        % So that f(E) = N*f_nodes
        % E = [nE nDims] , N = [nE nNodes]
        function N = evalAt(~,E)
            N = NaN(size(E,1),0) ;
        end
    end
    
% INTEGRATION QUADRATURE
    properties
        GaussIntegrationPoints = [] % [nGaussIntPts nDims]
        GaussIntegrationWeights = [] % [nGaussIntPts 1]
    end
    
%% CONSTRUCTOR / DESTRUCTOR
    methods
        function this = EmptyElement(varargin)
        % Constructor
        end

        function delete(this)
        % Destructor
        end
    end

end

