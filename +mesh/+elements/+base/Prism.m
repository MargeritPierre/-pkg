classdef Prism < pkg.mesh.elements.base.BaseElement
%PRISM 3D PRISM (2D simplex extruded along e3)

%% (MANDATORY) ELEMENT DEFINITION
% ELEMENT GEOMETRY
    properties (SetAccess = protected)
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates = [0 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1 ; 1 0 1 ; 0 1 1]
        % The list of faces [nFaces 4]
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
            N = [ (1-E(:,1)-E(:,2)).*(1-E(:,3)) E(:,1).*(1-E(:,3)) E(:,2).*(1-E(:,3)) ...
                  (1-E(:,1)-E(:,2)).*E(:,3) E(:,1).*E(:,3) E(:,2).*E(:,3) ] ;
        end
        
        function DER = evalDerivativeAt(this,E,ORD,~)
        % Evaluate the partial derivatives of order ORD = [o1 o2 ...] (:[1 nDims]) of 
        % shape functions at given local coordinates E:[nRows nDims]
        % exemple: df(E)/(dx²dy) = DER(E,[2 1 0]) ;
            if nargin<3 ; ORD = [1 0 0] ; end
            if numel(ORD)~=3 ; error('Wrong derivation order argument (must be [1 nDims])') ; end
            if all(ORD==0) % No Derivative
                DER = this.evalAt(E) ;
            elseif any(ORD>1) % High-order derivatives vanish
                DER = zeros(size(E,1),this.nNodes) ;
            else % Valid derivatives
                if isequal(ORD,[1 0 0])
                    DER = [ (1-E(:,3)).*[-1 1 0] E(:,3).*[-1 1 0]  ] ;
                elseif isequal(ORD,[0 1 0])
                    DER = [ (1-E(:,3)).*[-1 0 1] E(:,3).*[-1 0 1]  ] ;
                elseif isequal(ORD,[0 0 1])
                    DER = [ -(1-E(:,1)-E(:,2)) -E(:,1) -E(:,2) (1-E(:,1)-E(:,2)) E(:,1) E(:,2) ] ;
                elseif isequal(ORD,[1 1 0])
                    DER = zeros(size(E,1),this.nNodes) ;
                elseif isequal(ORD,[1 0 1])
                    DER = ones(size(E,1),1).*[1 -1 0 -1 1 0] ;
                elseif isequal(ORD,[0 1 1])
                    DER = ones(size(E,1),1).*[1 0 -1 -1 0 1] ;
                elseif isequal(ORD,[1 1 1])
                    DER = zeros(size(E,1),this.nNodes) ;
                end
            end
        end
    end
    
% INTEGRATION QUADRATURE
    properties
        GaussIntegrationPoints = [1/3 1/3 1.5] % [nGaussIntPts nDims]
        GaussIntegrationWeights = 1/2 % [nGaussIntPts 1]
    end
    
%% CONSTRUCTOR / DESTRUCTOR
    methods
        function this = Prism(varargin)
        % Constructor
            this.Edges = pkg.mesh.elements.ElementTable('Types',pkg.mesh.elements.base.Bar...
                                                            ,'Indices',[1 1 2 ; ...
                                                                        1 2 3 ; ...
                                                                        1 3 1 ; ...
                                                                        1 4 5 ; ...
                                                                        1 5 6 ; ...
                                                                        1 6 4 ; ...
                                                                        1 1 4 ; ...
                                                                        1 2 5 ; ...
                                                                        1 3 6]) ;
            this.Faces = pkg.mesh.elements.ElementTable('Types',[pkg.mesh.elements.base.Triangle ; pkg.mesh.elements.base.Quadrangle]...
                                                            ,'Indices',[1 1 3 2 0 ; ...
                                                                        1 4 5 6 0 ; ...
                                                                        2 1 2 5 4 ; ...
                                                                        2 2 3 6 5 ; ...
                                                                        2 3 1 4 6 ]) ;
        end

        function delete(this)
        % Destructor
        end
    end

end

