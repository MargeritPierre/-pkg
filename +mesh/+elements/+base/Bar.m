classdef Bar < pkg.mesh.elements.base.BaseElement
%BAR 1D element with two nodes

%% (MANDATORY) ELEMENT DEFINITION
% ELEMENT GEOMETRY
    properties (SetAccess = protected)
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates = [0 ; 1]
        % No faces
        Faces = pkg.mesh.elements.ElementTable
        % The only edge is itself (see constructor)
        Edges 
    end
    
% SHAPE FUNCTIONS AND DERIVATIVES
    methods
        % Evaluate the shape functions matrix N at local coordinates E
        % So that f(E) = N*f_nodes
        % E = [nE nDims==1] , N = [nE nNodes==2]
        function N = evalAt(~,E)
            N = [1-E(:) E(:)] ;
        end
        
        function DER = evalDerivativeAt(this,E,ORD,~)
        % Evaluate the partial derivatives of order ORD = [o1 o2 ...] (:[1 nDims]) of 
        % shape functions at given local coordinates E:[nRows nDims]
        % exemple: df(E)/(dx²dy) = DER(E,[2 1 0]) ;
            if nargin<3 ; ORD = 1; end
            if numel(ORD)~=1 ; error('Wrong derivation order argument (must be scalar)') ; end
            if ORD==0 % No Derivative
                DER = this.evalAt(E) ;
            elseif ORD>1 % High-order derivatives vanish
                DER = zeros(size(E,1),2) ;
            else % Valid derivative
                DER = ones(size(E,1),1).*[-1 1] ;
            end
        end
    end
    
% INTEGRATION QUADRATURE
    properties
        GaussIntegrationPoints = 0.5 % [nGaussIntPts nDims]
        GaussIntegrationWeights = 1 % [nGaussIntPts 1]
    end
    
%% CONSTRUCTOR / DESTRUCTOR
    methods
        function this = Bar(varargin)
        % Constructor
            this.Edges = pkg.mesh.elements.ElementTable('Types',this,'Indices',[1 1 2]) ;
        end

        function delete(this)
        % Destructor
        end
    end
    
%% GEOMETRY
    methods
        function [elems,idx] = slice(this,nodeBool)
        % Slice the element given a signed logical value of a levelset on
        % each node: -1 (inside), 0 (on) or 1 (outside)
        % input: this: element; nodeBool [nElems this.nNodes]
        % output:
        %   - an element table ELEMS containing the sliced elements. 
        %       /!\ the node indices in the table are complex uint32 ! :
        %       - real indices denote nodes of the reference element
        %       - imaginary indices denote edges of the reference element
        %   - a list IDX of size [ELEMS.nElems 1] where IDX(i) contains the
        %   index of the input element
            ON = nodeBool==0 ;
            elems = pkg.mesh.elements.ElementTable ;
            idx = [] ;
        % If the two nodes are ON the levelset, then return the same element
            bool = all(ON,2) ;
            elems = [elems ; ...
                        pkg.mesh.elements.ElementTable(...
                                    'Types',this ...
                                    ,'Indices',repmat([1 1 2],[sum(bool) 1]) ...
                        ) ...
                    ] ;
            idx = [idx(:) ; find(bool)] ;
        % If only one of the two nodes is ON the levelset, then return a
        % node on the edge
            bool = sum(ON,2)==1 ;
            elems = [elems ; ...
                        pkg.mesh.elements.ElementTable(...
                                    'Types',pkg.mesh.elements.base.Node ...
                                    ,'Indices',repmat([1 1i],[sum(bool) 1]) ...
                        ) ...
                    ] ;
            idx = [idx(:) ; find(bool)] ;
        end
    end

end

