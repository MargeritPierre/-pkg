classdef TableObject < handle
%TABLEOBJECT Interface for an object containing a table of elements
% Used in pkg.geometry.mesh.Mesh, pkg.geometry.mesh.FunctionSpace, etc.
    
properties (SetObservable)
% The table of elements
    Elems = pkg.geometry.mesh.elements.ElementTable
end
properties (SetAccess = protected)
% The table of faces (semi-dependent)
    Faces pkg.geometry.mesh.elements.ElementTable = pkg.geometry.mesh.elements.ElementTable
% The table of edges (semi-dependent)
    Edges pkg.geometry.mesh.elements.ElementTable = pkg.geometry.mesh.elements.ElementTable
end
properties (Hidden)
    % Elem to features sparse matrices
    ElemFaces
    ElemEdges
end


%% CONSTRUCTOR/DESTRUCTOR
methods
    function this = TableObject(table)
    % Constructor
        if nargin==0 ; return ; end
        this.Elems = table ;
    end
    
    function delete(this)
    % Destructor
    end
end


%% COUNT FUNCTIONS
methods
    function N = nElems(this) ; N = this.Elems.nElems ; end
    function N = nEdges(this) ; N = this.Edges.nElems ; end
    function N = nFaces(this) ; N = this.Faces.nElems ; end
    function N = nMaxNodesByElem(this) ; N = this.Elems.nMaxNodesByElem ; end
    function N = nNodesInElem(this) ; N = this.Elems.nNodes ; end
end
    

%% ELEMENT SET/GET INTERFACE
methods
    function set.Elems(this,elems)
    % Modify the table of elements (create a new table each time)
    % Depending on the provided data...
        switch class(elems)
            case 'pkg.geometry.mesh.elements.ElementTable'
            otherwise
                if ~ismatrix(elems) ; error('wrong format: must be [nElems nMaxNodesByElems]') ; end
                elems = uint32(elems) ;
                if ~isempty(this.Elems.Types)
                    elems = pkg.geometry.mesh.elements.ElementTable('Types',this.Elems.Types,'Indices',elems) ;
                else
                    elems = pkg.geometry.mesh.elements.ElementTable('Types',this.Elems.Types,'NodeIdx',elems) ;
                end
        end
    % Erase mesh feature lists | 
    % will be computed the first time its needed
    % (see get.Faces/Edges below)
        this.Faces = pkg.geometry.mesh.elements.ElementTable ; 
        this.ElemFaces = [] ;
        this.Edges = pkg.geometry.mesh.elements.ElementTable ; 
        this.ElemEdges = [] ;
    % Set
        this.Elems = elems ;
    end
    
    % Compute the faces & edges only when needed
    function edges = get.Edges(this)
        computeFeatures(this,'Edges',this.Edges) ;
        edges = this.Edges ;
    end
    function elmedg = get.ElemEdges(this)
        computeFeatures(this,'Edges',this.Edges) ;
        elmedg = this.ElemEdges ;
    end
    function faces = get.Faces(this)
        computeFeatures(this,'Faces',this.Faces) ;
        faces = this.Faces ;
    end
    function elmfac = get.ElemFaces(this)
        computeFeatures(this,'Faces',this.Faces) ;
        elmfac = this.ElemFaces ;
    end
    function computeFeatures(this,feat,data)
        if ~isempty(this.Elems) && isempty(data) 
            [this.(feat),this.(['Elem' feat])] = this.Elems.getTableOfUnique(feat) ; 
        end
    end
end


%% CONNECTIVITY MATRICES
methods
    function M = elem2node(this,varargin)
    % Elements to Nodes [nNodes nElems]: nod(:) = M*elmt(:)
        M = sparse(this.Elems,varargin{:}) ;
    end
    function M = face2node(this,varargin)
    % Faces to Nodes [nNodes nFaces]: nod(:) = M*fa(:)
        M = sparse(this.Faces,varargin{:}) ;
    end
    function M = edge2node(this,varargin)
    % Faces to Nodes [nNodes nFaces]: nod(:) = M*fa(:)
        M = sparse(this.Edges,varargin{:}) ;
    end
    function M = elem2face(this,varargin)
    % Elements to Edges [nFaces nElems]: edg(:) = M*elmt(:)
        M = logical(this.ElemFaces) ;
    end
    function M = elem2edge(this,varargin)
    % Elements to Edges [nEdges nElems]: edg(:) = M*elmt(:)
        M = logical(this.ElemEdges) ;
    end
    function M = face2edge(this,varargin)
    % Faces to Edges [nEdges nFaces]: edg(:) = M*face(:)
        M = this.Faces.contains(this.Edges) ;
    end
end


%% AUTOMATIC ELEMENT ASSIGNMENT
methods
    function assignElements(this,types)
    % Automatically assign elements types in the table
        if nargin<2 ; types = '2D' ; end
    % Process types options
        if ischar(types)
            switch types
                case '1D'
                    geo = {'Node' 'Bar' } ;
                case '2D'
                    geo = {'Node' 'Bar' 'Triangle' 'Quadrangle'} ;
                case '3D'
                    geo = {'Node' 'Bar' 'Tetrahedron' 'Hexahedron' 'Prism' 'Pyramid'} ;
                otherwise
                    error('Unsupported option for element types.') ;
            end
            types = pkg.geometry.mesh.elements.AbstractElement.empty ;
            for tt = 1:numel(geo) ; types(end+1) = eval(['pkg.geometry.mesh.elements.base.' geo{tt}]) ; end
        end
    % Count nodes in each existing element
        nNodesInElems = this.Elems.nNodes ;
    % Check that types match node numbers
        valid = arrayfun(@(t)ismember(t.nNodes,nNodesInElems),types) ;
    % Assign
        this.Elems.Types = types(valid) ;
    end
end



%% SPECIAL FEATURES
methods
    
    function ends = endNodes(this)
    % Return an array of logical, true if the node is an end node
    % A END node belongs to one or less elements
        ends = sum(this.elem2node,2)<=1 ;
    end

    function outFace = outerFaces(this)
    % Return an array of logical, true if the face is one an the outer surface
    % The face is on an outer surface if it is linked to one element at most
        outFace = sum(this.elem2face,2)<=1 ;
    end

    function outNod = outerNodes(this)
    % Return an array of logical, true if the node is on an the outer surface
    % The node is on an outer surface if it is linked to an outer edge or
    % is an end node
        outNod = logical(this.edge2node*outerEdges(this)) ;
        outNod = outNod | this.endNodes ;
    end

    function outEdg = outerEdges(this)
    % Return an array of logical, true if the edge is on the outer surface
    % The edge is on the surface if it is linked to an outer face
    % or is a boundary edge
        outEdg = logical(this.face2edge*outerFaces(this)) ;
        outEdg = outEdg | this.boundaryEdges ;
    end

    function bndEdg = boundaryEdges(this)
    % Return an array of logical, true if the edge is on the boundary
    % The edge is on the boundary if it is linked to at most one element
        bndEdg = sum(this.elem2edge,2)<=1 ;
    end

    function bndNod = boundaryNodes(this)
    % Return an array of logical, true if the node is on the boundary
    % The node is on the boundary if it is linked to a boundary edge
        bndEdg = boundaryEdges(this) ;
        bndNod = logical(this.edge2node*bndEdg) ;
    end

    function [crv,asCell,edgBndCrv] = boundaryCurves(this)
    % Return an array of node indices representing the sorted mesh boundaries
    % crv contains the nodes indices with the different curves
    % separated by NaNs
    % asCells is the same but the node indices are in different cells
    % crvEdg [1 this.nEdges] uint32 conatins the index of the bnd crv
    % associated to each edge (==0 if interior edge) 
        if this.nElems==0 ; crv = {} ; return ; end
        % List of boundary edges
            bndEdg = find(boundaryEdges(this)) ;
            nn = this.Edges.NodeIdx(bndEdg,:) ; % boundary edge nodes
        % Scan the edges to build the curves
            edgBndCrv = zeros(1,this.nEdges,'uint32') ;
            crv{1} = nn(1,:) ;
            nn(1,:) = [] ;
            edgBndCrv(bndEdg(1)) = 1 ; 
            bndEdg(1) = [] ;
            while ~isempty(nn)
                % Find the next edge and associated end node
                    [nextEdg,indNode] = find(nn==crv{end}(end),1,'first') ;
                % If no next edge found, create a new curve
                    if isempty(nextEdg)
                        crv{end+1} = [] ;
                        nextEdg = 1 ;
                        indNode = 0 ;
                    end
                % Remaining edge nodes that have to be added
                    if indNode==0 ; addNodes = nn(nextEdg,1:end) ;
                    elseif indNode==1 ; addNodes = nn(nextEdg,2:end) ;
                    else ; addNodes = nn(nextEdg,indNode-1:-1:1) ; 
                    end
                    addNodes = addNodes(addNodes~=0) ;
                % Add the nodes to the current curve
                    crv{end} = [crv{end} addNodes] ;
                    nn(nextEdg,:) = [] ;
                % Set the edge boundary curve number
                    edgBndCrv(bndEdg(nextEdg)) = numel(crv) ; 
                    bndEdg(nextEdg) = [] ;
            end
        % Concatenate curves
            crv = reshape(crv,1,[]) ;
            asCell = crv ; % Backup individual curves
            crv(end+1,:) = {NaN} ;
            crv = cat(2,crv{:}) ;
            crv = crv(1:end-1) ;
    end

end

end

