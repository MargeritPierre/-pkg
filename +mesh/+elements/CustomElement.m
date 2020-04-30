classdef CustomElement < pkg.mesh.elements.AbstractElement
%CUSTOMELEMENT Template for a custom element class
% Takes the example of a 9-nodes P2 quad
%
% New Element class creation: 
%   - Copy-paste this file
%   - Rename it with the newclass name "myElementName"
%   - Find "CustomElement" / Replace by "myElementName"
%   - Modify the properties and methods to fit you needs

%% (MANDATORY) ELEMENT DEFINITION
% ELEMENT GEOMETRY
    properties (SetAccess = protected)
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates = [-1 -1 ; 0 -1 ; 1 -1 ; 1 0 ; 1 1 ; 0 1 ; -1 1 ; -1 0 ; 0 0]
        % The list of faces [nFaces nMaxNodesByFace]
        %Faces = pkg.mesh.elements.ElementTable([1 2 3 4 5 6 7 8]) % nineth node is slave
        Faces = pkg.mesh.elements.ElementTable([1 2 9 8 ; 2 3 4 9 ; 9 4 5 6 ; 8 9 6 7]) % nineth node is included
        % The list of edges [nEdges 2] 
        Edges =  pkg.mesh.elements.ElementTable([1 2 ; 2 3 ; 3 4 ; 4 5 ; 5 6 ; 6 7 ; 7 8 ; 8 1])
    end
    
% SHAPE FUNCTIONS AND DERIVATIVES
    methods
        % Evaluate the shape functions matrix N at local coordinates E
        % So that f(E) = N*f_nodes
        % E = [nE nDims] , N = [nE nNodes]
        function N = evalAt(~,E)
        % See p.21 of iut.univ-lemans.fr/ydlogi/cours/mef_elas_2d.pdf
        % L = [L1(e1) L1(e2) L2(e1) L2(e2) L3(e1) L3(e2)]
            L = [E.*(E-1)/2 1-E.^2 E.*(E+1)/2] ;
            N = [L(:,1).*L(:,2) , ...
                 L(:,3).*L(:,2) , ...
                 L(:,5).*L(:,2) , ...
                 L(:,5).*L(:,4) , ...
                 L(:,5).*L(:,6) , ...
                 L(:,3).*L(:,6) , ...
                 L(:,1).*L(:,6) , ...
                 L(:,1).*L(:,4) , ...
                 L(:,3).*L(:,4) ] ;
        end
        
        function DER = evalDerivativeAt(this,E,ORD,~)
        % Evaluate the partial derivatives of order ORD = [o1 o2 ...] (:[1 nDims]) of 
        % shape functions at given local coordinates E:[nRows nDims]
        % When possible, it is preferred to override this method to provide
        % analytical derivatives !
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
        GaussIntegrationPoints = [-1 -1 ; 1 -1 ; 1 1 ; -1 1] % [nGaussIntPts nDims]
        GaussIntegrationWeights = 1/4 * [1 ; 1 ; 1 ; 1 ] % [nGaussIntPts 1]
    end
    
%% CONSTRUCTOR / DESTRUCTOR
    methods
        function this = CustomElement(varargin)
        % Constructor
        end

        function delete(this)
        % Destructor
        end
    end

end

