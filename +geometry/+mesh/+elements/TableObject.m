classdef TableObject < matlab.mixin.Copyable
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
    function M = feat2node(this,M)
    % Adjust the size of the connectivity matrix to match the number of nodes
        if ~ismethod(this,'nNodes') ; return ; end
        M(end+1:this.nNodes,:) = 0 ; 
    end
    function M = elem2node(this,varargin)
    % Elements to Nodes [nNodes nElems]: nod(:) = M*elmt(:)
        M = this.feat2node(sparse(this.Elems,varargin{:})) ;
    end
    function M = face2node(this,varargin)
    % Faces to Nodes [nNodes nFaces]: nod(:) = M*fa(:)
        M = this.feat2node(sparse(this.Faces,varargin{:})) ;
    end
    function M = edge2node(this,varargin)
    % Faces to Nodes [nNodes nFaces]: nod(:) = M*fa(:)
        M = this.feat2node(sparse(this.Edges,varargin{:})) ;
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
% END FEATURES
    function [no,el,ed,fa] = endFeatures(this)
    % Return end features
        el2no = this.elem2node ;
        ed2no = this.edge2node ;
    % End nodes
        no = (sum(el2no,2)==0) ... belongs to no element
             | (sum(el2no,2)==1 & sum(ed2no,2)>0) ... belong to an edge of ONE element
             ;
     % Other features
        if nargout>=2 ; el = logical(el2no'*no) ; end
        if nargout>=3 ; ed = logical(ed2no'*no) ; end
        if nargout>=4 ; fa = logical(this.face2edge'*ed) ; end
    end
    
    function no = endNodes(this)
    % Return an array of logical, true if the node is an end node
    % A END node belongs to one or less elements
        no = endFeatures(this) ;
    end
    
    function el = endElems(this)
    % Return an array of logical, true if the elem is linked to an end node
        [~,el] = endFeatures(this) ;
    end
    
    function ed = endEdges(this)
    % Return an array of logical, true if the edge is linked to an end node
        [~,~,ed] = endFeatures(this) ;
    end
    
    function fa = endFaces(this)
    % Return an array of logical, true if the face is linked to an end node
        [~,~,~,fa] = endFeatures(this) ;
    end
    
% OUTER FEATURES
    function [fa,el,ed,no] = outerFeatures(this)
    % Return outer features: attached to a outer face
        e2f = this.elem2face ;
        fa = sum(e2f,2)<=1 ;
        if nargout>=2 ; el = logical(e2f'*fa) ; end
        if nargout>=3 ; ed = logical(this.face2edge*fa) | this.boundaryEdges ; end
        if nargout>=4 ; no = logical(this.edge2node*ed) | this.endNodes ; end
    end
    
    function fa = outerFaces(this)
    % Return an array of logical, true if the face belongs to the outer surface
    % The face is on an outer surface if it is linked to one element at most
        fa = this.outerFeatures ;
    end
    
    function el = outerElems(this)
    % Return an array of logical, true if the elem is on the outer surface
    % The elem is on an outer surface if it is linked to one outer face
        [~,el] = this.outerFeatures ;
    end
    
    function ed = outerEdges(this)
    % Return an array of logical, true if the edge is on the outer surface
    % The edge is on an outer surface if it is linked to one outer face or
    % is a boundary edge
        [~,~,ed] = this.outerFeatures ;
    end

    function no = outerNodes(this)
    % Return an array of logical, true if the node is on an the outer surface
    % The node is on an outer surface if it is linked to an outer edge or
    % is an end node
        [~,~,~,no] = this.outerFeatures ;
    end

% BOUNDARY FEATURES
    function [ed,el,no,fa] = boundaryFeatures(this)
    % Return features associated to boundaries
        ele2edg = this.elem2edge ;
        ed = sum(ele2edg,2)<=1 ;
        if nargout>=2 ; el = logical(ele2edg'*ed) ; end
        if nargout>=3 ; no = logical(this.edge2node*ed) | this.endNodes ; end
        if nargout>=4 ; fa = logical(this.face2edge'*ed) ; end
    end
    
    function ed = boundaryEdges(this)
    % Return an array of logical, true if the edge is on the boundary
    % The edge is on the boundary if it is linked to at most one element
        ed = boundaryFeatures(this) ;
    end

    function el = boundaryElems(this)
    % Return an array of logical, true if the element is on the boundary
    % The element is on the boundary if it is linked to a boundary edge
        [~,el] = boundaryFeatures(this) ;
    end

    function no = boundaryNodes(this)
    % Return an array of logical, true if the node is on the boundary
    % The node is on the boundary if it is linked to a boundary edge or is
    % an end node
        [~,~,no] = boundaryFeatures(this) ;
    end

    function fa = boundaryFaces(this)
    % Return an array of logical, true if the face is on the boundary
    % The face is on the boundary if it is linked to a boundary edge
        [~,~,~,fa] = boundaryFeatures(this) ;
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

