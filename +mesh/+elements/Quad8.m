classdef Quad8 < pkg.mesh.elements.AbstractElement
%QUAD8 element following the ABAQUS numerotation
% see http://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node43.html

%% (MANDATORY) ELEMENT DEFINITION
% ELEMENT GEOMETRY
    properties (SetAccess = protected)
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates = [0 0 ; 1 0 ; 1 1 ; 0 1 ; 0.5 0 ; 1 0.5 ; 0.5 1 ; 0 0.5]
        % The list of faces (see constructor below)
        Faces % (==this)
        % The list of edges (see constructor below)
        Edges  
    end
    
% SHAPE FUNCTIONS AND DERIVATIVES
    methods
        % Evaluate the shape functions matrix N at local coordinates E
        % So that f(E) = N*f_nodes
        % E = [nE nDims] , N = [nE nNodes]
        function N = evalAt(~,E)
        % See p.20 of iut.univ-lemans.fr/ydlogi/cours/mef_elas_2d.pdf
        % Map the usual [-1;1]x[-1;1] domain of this element to the current [0;1]x[0;1]
            E = E*2-1 ;
        % Shape functions
            N = [-(1-E(:,1)).*(1-E(:,2)).*(1+E(:,1)+E(:,2)) , ...
                 -(1+E(:,1)).*(1-E(:,2)).*(1-E(:,1)+E(:,2)) , ...
                 -(1+E(:,1)).*(1+E(:,2)).*(1-E(:,1)-E(:,2)) , ...
                 -(1-E(:,1)).*(1+E(:,2)).*(1+E(:,1)-E(:,2)) , ...
                 2.*(1-E(:,1).^2).*(1-E(:,2)) , ...
                 2.*(1+E(:,1)).*(1-E(:,2).^2) , ...
                 2.*(1-E(:,1).^2).*(1+E(:,2)) , ...
                 2.*(1-E(:,1)).*(1-E(:,2).^2) ...
                 ] ;
        end
        
%         function DER = evalDerivativeAt(this,E,ORD,~)
%         % Evaluate the partial derivatives of order ORD = [o1 o2 ...] (:[1 nDims]) of 
%         % shape functions at given local coordinates E:[nRows nDims]
%         % When possible, it is preferred to override this method to provide
%         % analytical derivatives !
%         % exemple: df(E)/(dx²dy) = DER(E,[2 1 0]) ;
%             if nargin<3 ; ORD = [1 zeros(1,this.nDims-1)] ; end
%             if numel(ORD)~=this.nDims ; error('Wrong derivation order argument (must be [1 nDims])') ; end
%             if all(ORD==0) % No Derivative
%                 DER = this.evalAt(E) ;
%             elseif any(ORD>2) % High-order derivatives vanish
%                 DER = zeros(size(E,1),this.nNodes) ;
%             else % Valid derivatives
%                 % See p.21 of iut.univ-lemans.fr/ydlogi/cours/mef_elas_2d.pdf
%                 % L = [L1(e1) L1(e2) L2(e1) L2(e2) L3(e1) L3(e2)]
%                 DL = [E.*(E-1)/2 1-E.^2 E.*(E+1)/2] ; % zeroth order
%                 DL(:,:,end+1) = [E-1/2 -2*E E+1/2] ; % first order
%                 if any(ORD==2) ; DL(:,:,end+1) = repmat([1 1 -2 -2 1 1],[size(E,1) 1]) ; end % second order
%                 % Apply the derivative
%                 DER = DL(:,[1 3 5 5 5 3 1 1 3],ORD(1)+1).*DL(:,[2 2 2 4 6 6 6 4 4],ORD(2)+1) ;
%             end
%         end
    end
    
% INTEGRATION QUADRATURE
    properties
        GaussIntegrationPoints = [-1 -1 ; 1 -1 ; 1 1 ; -1 1] % [nGaussIntPts nDims]
        GaussIntegrationWeights = 1/4 * [1 ; 1 ; 1 ; 1 ] % [nGaussIntPts 1]
    end
    
%% CONSTRUCTOR / DESTRUCTOR
    methods
        function this = Quad8(varargin)
        % Constructor
            this = this@pkg.mesh.elements.AbstractElement(varargin{:}) ;
            this.Faces = pkg.mesh.elements.ElementTable('Types',this,'NodeIdx',1:8) ;
            this.Edges = pkg.mesh.elements.ElementTable('Types',pkg.mesh.elements.LagrangeElement('1D',2),'NodeIdx',[1 5 2 ; 2 6 3 ; 3 7 4 ; 4 8 1]) ;
        end

        function delete(this)
        % Destructor
        end
    end

end

