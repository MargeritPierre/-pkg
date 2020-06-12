classdef Quadrangle < pkg.mesh.elements.base.BaseElement
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
        % exemple: df(E)/(dx²dy) = DER(E,[2 1 0]) ;
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
    
%% GEOMETRY
    methods
        function [elems,idx] = slice(this,nodeBool)
        % Slice the element given a signed logical value of a levelset on
        % each node: -1 (inside), 0 (on) or 1 (outside)
        % input: this: element; nodeBool [nElems this.nNodes]
        % output:
        %   - an element table ELEMS containing the sliced elements. 
        %       /!\ the node indices in the table are COMPLEX uint32 ! :
        %       - indices<=this.nNodes denote nodes of the reference element
        %       - indices>this.nNodes denote edges of the reference element
        %           (edgIdx = idx-this.nNodes)
        %   - a list IDX of size [ELEMS.nElems 1] where IDX(i) contains the
        %   index of the input element
            ON = nodeBool==0 ;
            elems = pkg.mesh.elements.ElementTable ;
            idx = [] ;
        % If all nodes are ON the levelset, then return the same element
            bool = all(ON,2) ;
            if any(bool)
                elems = [elems ; ...
                            pkg.mesh.elements.ElementTable(...
                                        'Types',pkg.mesh.elements.base.Quadrangle ...
                                        ,'Indices',repmat([1 1 2 3 4],[sum(bool) 1]) ...
                            ) ...
                        ] ;
                idx = [idx(:) ; find(bool)] ;
            end
        % Only one node ON the levelset: add a node
            bool = sum(ON,2)==1 ;
            if any(bool)
                [~,aloneNode] = find(ON(bool,:)) ;
                elems = [elems ; ...
                            pkg.mesh.elements.ElementTable(...
                                        'Types',pkg.mesh.elements.base.Node ...
                                        ,'Indices',aloneNode(:).*[0 1] + [1 0] ...
                            ) ...
                        ] ;
                idx = [idx(:) ; find(bool)] ;
            end
        % Exactly two nodes ON the levelset: add a bar joining the two
            bool = sum(ON,2)==2 ;
            if any(bool)
                [~,ind] = sort(ON(bool,:),2,'descend') ;
                elems = [elems ; ...
                            pkg.mesh.elements.ElementTable(...
                                        'Types',pkg.mesh.elements.base.Bar ...
                                        ,'Indices',[ones(size(ind,1),1) ind(:,1:2)] ...
                            ) ...
                        ] ;
                idx = [idx(:) ; find(bool)] ;
            end
        % Exactly three nodes ON the levelset: add a triangle
            bool = sum(ON,2)==3 ;
            if any(bool)
                [~,ind] = sort(ON(bool,:),2,'descend') ;
                elems = [elems ; ...
                            pkg.mesh.elements.ElementTable(...
                                        'Types',pkg.mesh.elements.base.Triangle ...
                                        ,'Indices',[ones(size(ind,1),1) ind(1:3)] ...
                            ) ...
                        ] ;
                idx = [idx(:) ; find(bool)] ;
            end
        % EDGE CROSSING
            % One node is out/in, all other nodes are in/out
                config = [1 -1 -1 -1] ;
                for cc = 0:3
                    bTest = circshift(config,cc,2) ;
                    nodEdges = circshift(1:4,-cc,2) ;
                    for ss = [-1 1]
                        bool = all(nodeBool==bTest.*ss,2) ;
                        if any(bool)
                            elems = [elems ; ...
                                        pkg.mesh.elements.ElementTable(...
                                                    'Types',pkg.mesh.elements.base.Bar ...
                                                    ,'Indices',repmat([1 this.nNodes+[nodEdges([end,1])]],[sum(bool) 1]) ...
                                        ) ...
                                    ] ;
                            idx = [idx(:) ; find(bool)] ;
                        end
                    end
                end
            % Two adjacent nodes are out/in, the two other nodes are in/out
                config = [1 1 -1 -1] ;
                for cc = 0:3
                    bTest = circshift(config,cc,2) ;
                    nodEdges = circshift(1:4,-cc,2) ;
                    for ss = [-1 1]
                        bool = all(nodeBool==bTest.*ss,2) ;
                        if any(bool)
                            elems = [elems ; ...
                                        pkg.mesh.elements.ElementTable(...
                                                    'Types',pkg.mesh.elements.base.Bar ...
                                                    ,'Indices',repmat([1 this.nNodes+[nodEdges([end,2])]],[sum(bool) 1]) ...
                                        ) ...
                                    ] ;
                            idx = [idx(:) ; find(bool)] ;
                        end
                    end
                end
        end
        
%         function S = sliceCases(this)
%         % Return all the cases that have to be tested while slicing this
%         % element
%         % S is an array of structures with fields:
%         %   - Test: signed boolean configuration [1 this.nNodes]
%         %           (-1 inside, 0 on, 1 outside, NaN any cases)
%         %   - IN: elements inside (pkg.mesh.elements.ElementTable)
%         %   - ON: elements on the slice (pkg.mesh.elements.ElementTable)
%         %   - OUT: elements outside (pkg.mesh.elements.ElementTable)
%         % /!\ The parent function will already test opposite signs, do not bother with this !
%         % in the ElementTables, indices>this.nNodes denote edges indices
%         % (used when the slice cross an edge)
%             import pkg.mesh.elements.ElementTable
%             import pkg.mesh.elements.base.*
%             S = repmat(struct('Test',[],'IN',[],'ON',[],'OUT',[]),[0 0]) ; % empty structure
%         % One node on the slice, all other nodes outside
%             for nn = 1:4
%                 S(end+1).Test = circshift([0 1 1 1],nn-1,2) ;
%                 S(end).ON = ElementTable('Types',Node,'Indices',[1 nn]) ;
%                 S(end).OUT = ElementTable('Types',Quadrangle,'Indices',[1 1 2 3 4]) ;
%             end
%         % Two adjacent nodes on the slice, the two other point outside
%             for nn = 1:4
%                 S(end+1).Test = circshift([0 0 1 1],nn-1,2) ;
%                 S(end).ON = ElementTable('Types',Bar,'Indices',[1 find(S(end).Test==0)]) ;
%                 S(end).OUT = ElementTable('Types',Quadrangle,'Indices',[1 1 2 3 4]) ;
%             end
%         % Two adjacent nodes on the slice, other nodes on opposite sides
%             for nn = 1:4
%                 S(end+1).Test = circshift([0 0 1 -1],nn-1,2) ;
%                 edges = circshift(1:4,1-nn,2) ;
%                 idx = find(S(end).Test==0).*[1 1;1 0 ;0 1] + [0 0 ; 0 1 ; 1 0]*(edges(3)+this.nNodes) ;
%                 S(end).ON = ElementTable('Types',Bar,'Indices',[[1;1;1] idx]) ;
%                 S(end).OUT = ElementTable('Types',Quadrangle,'Indices',[1 1 2 3 4]) ;
%             end
%         end
    end

end

