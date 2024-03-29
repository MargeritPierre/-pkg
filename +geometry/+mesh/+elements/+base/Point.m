classdef Point < pkg.geometry.mesh.elements.base.BaseElement
%Point 0D element representing a lone node

%% (MANDATORY) ELEMENT DEFINITION
% ELEMENT GEOMETRY
    properties (SetAccess = protected)
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates = zeros(1,0)
        % No faces
        Faces = pkg.geometry.mesh.elements.ElementTable
        % The only edge is itself (see constructor)
        Edges = pkg.geometry.mesh.elements.ElementTable
    end
    
% SHAPE FUNCTIONS AND DERIVATIVES
    methods
        % Evaluate the shape functions matrix N at local coordinates E
        % So that f(E) = N*f_nodes
        % E = [nE nDims==1] , N = [nE nNodes==2]
        function N = evalAt(~,E)
            N = ones(size(E,1),1) ;
        end
        
        function DER = evalDerivativeAt(this,E,ORD,~)
        % Evaluate the partial derivatives of order ORD = [o1 o2 ...] (:[1 nDims]) of 
        % shape functions at given local coordinates E:[nRows nDims]
        % exemple: df(E)/(dx�dy) = DER(E,[2 1 0]) ;
            if nargin<3 ; ORD = 1; end
            if numel(ORD)~=1 ; error('Wrong derivation order argument (must be scalar)') ; end
            if ORD==0 % No Derivative
                DER = this.evalAt(E) ;
            else
                DER = zeros(size(E,1),1) ;
            end
        end
    end
    
% INTEGRATION QUADRATURE
    properties
        GaussIntegrationPoints = [] % [nGaussIntPts nDims]
        GaussIntegrationWeights = [] % [nGaussIntPts 1]
    end
    
%% CONSTRUCTOR / DESTRUCTOR
    methods
        function this = Point(varargin)
        % Constructor
        end

        function delete(this)
        % Destructor
        end
    end
    
%% GEOMETRY
    methods
%         function [elems,idx] = slice(this,nodeBool)
%         % Slice the element given a signed logical value of a levelset on
%         % each node: -1 (inside), 0 (on) or 1 (outside)
%         % input: this: element; nodeBool [nElems this.nNodes]
%         % output:
%         %   - an element table ELEMS containing the sliced elements. 
%         %       /!\ the node indices in the table are complex uint32 ! :
%         %       - real indices denote nodes of the reference element
%         %       - imaginary indices denote edges of the reference element
%         %   - a list IDX of size [ELEMS.nElems 1] where IDX(i) contains the
%         %   index of the input element
%             ON = nodeBool==0 ;
%             elems = pkg.geometry.mesh.elements.ElementTable ;
%             idx = [] ;
%         % If the two nodes are ON the levelset, then return the same element
%             bool = all(ON,2) ;
%             elems = [elems ; ...
%                         pkg.geometry.mesh.elements.ElementTable(...
%                                     'Types',this ...
%                                     ,'Indices',repmat([1 1 2],[sum(bool) 1]) ...
%                         ) ...
%                     ] ;
%             idx = [idx(:) ; find(bool)] ;
%         % If only one of the two nodes is ON the levelset, then return a
%         % node on the edge
%             bool = sum(ON,2)==1 ;
%             elems = [elems ; ...
%                         pkg.geometry.mesh.elements.ElementTable(...
%                                     'Types',pkg.geometry.mesh.elements.base.Node ...
%                                     ,'Indices',repmat([1 1i],[sum(bool) 1]) ...
%                         ) ...
%                     ] ;
%             idx = [idx(:) ; find(bool)] ;
%         end
    end

end

