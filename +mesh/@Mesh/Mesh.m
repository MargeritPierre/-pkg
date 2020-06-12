classdef Mesh < matlab.mixin.Copyable
% The general mesh class


%% CONSTRUCTOR / DESTRUCTOR / COPY / SAVE / LOAD
methods
    function this = Mesh(varargin)
    % Class Constructor
    % Initialize the node position
        this.X = pkg.mesh.MeshFunction('Mesh',this) ;
    % Process inputs
        if nargin==1 % mesh from single input (see below)
            this = this.meshFromOneInput(varargin{:}) ; 
            return ; 
        end
        if mod(nargin,2) ; error('wrong number of arguments') ; end
    % Input arguments
        PropNames = varargin(1:2:end-1) ;
        PropValues = varargin(2:2:end) ;
    % Compatibility fixes
        [PropNames{ismember(PropNames,{'Nodes'})}] = deal('X') ;
    % Process input arguments
        for pp = 1:numel(PropNames)
            this.(PropNames{pp}) = PropValues{pp} ;
        end
    % Default element types
        if isempty(this.Elems.Types)
            this.assignElements('auto') ;
        end
    end
    
    function this = meshFromOneInput(this,input)
    % Build meshes from a simple input 
    % Input is the node coordinates
        if isnumeric(input)
            x = input ;
            switch ndims(x)
                case 2 % wireframe mesh
                    nodeIdx = [1:size(x,1)-1 ; 2:size(x,1)]' ;
                    types = pkg.mesh.elements.base.Bar ;
                case 3 % surface mesh
                    [nY,nX,~] = size(x) ;
                    p1 = (1:nY-1)'+ nY*(0:nX-2) ; p2 = p1+nY ; p3 = p2+1 ; p4 = p1+1 ;
                    nodeIdx = [p1(:) p2(:) p3(:) p4(:)] ;
                    types = pkg.mesh.elements.base.Quadrangle ;
                case 4 % volume mesh
                    error('Automatic volume mesh not implemented yet') ;
            end
            this.X.Values = reshape(x,[],size(x,ndims(x))) ;
            this.Elems = pkg.mesh.elements.ElementTable('Types',types,'NodeIdx',nodeIdx) ;
        else
    % Mesh from a specific object class
            switch class(input)
                case 'pkg.mesh.Mesh'
                    this = copy(this) ;
                otherwise
                    error(['Cannot create a mesh from a ' class(input)])
            end
        end
    end
end
methods (Access = protected)
    function mesh = copyElement(this)
    % Retrun a copy of the object
        mesh = copyElement@matlab.mixin.Copyable(this) ;
        mesh.X = copy(this.X) ;
        mesh.X.Mesh = mesh ;
    end
end


%% MESH NODES & ELEMENTS
properties (SetObservable = true)
    % Nodal coordinates (position in 1-2-3D space)
    X
    % Table of elements 
    Elems = pkg.mesh.elements.ElementTable
end
methods
    function set.Elems(this,elems)
    % Modify the mesh table of elements (create a new table each time)
    % Depending on the privided data...
        switch class(elems)
            case 'pkg.mesh.elements.ElementTable'
            otherwise
                if ~ismatrix(elems) ; error('wrong format: must be [nElems nMaxNodesByElems]') ; end
                elems = uint32(elems) ;
                if ~isempty(this.Elems.Types)
                    elems = pkg.mesh.elements.ElementTable('Types',this.Elems.Types,'Indices',elems) ;
                else
                    elems = pkg.mesh.elements.ElementTable('Types',this.Elems.Types,'NodeIdx',elems) ;
                end
        end
    % Clean the element list
        %elems = clean(elems) ; % on mesh cleanup only
    % Erase mesh feature lists | 
    % will be computed the first time its needed
    % (see get.Faces/Edges below)
        this.Faces = [] ; %elems.getTableOfUnique('Faces') ;
        this.Edges = [] ; %elems.getTableOfUnique('Edges') ;
    % Set
        this.Elems = elems ;
    end

    function set.X(this,values)
    % Set the mesh node coordinates
    % handle values or Mesh Functions
        switch class(values)
            case 'pkg.mesh.MeshFunction'
                values.Mesh = this ;
                this.X = values ;
            otherwise
                this.X.Values = values ;
        end
    end
    
    function assignElements(this,types)
    % Automatically assign elements types to the mesh
        if nargin<2 ; types = 'auto' ; end
    % Process types options
        if ischar(types)
            switch types
                case 'auto' % Preset list of base elements
                    if isPlanar(this) ; geo = {'Node' 'Bar' 'Triangle' 'Quadrangle'} ;
                    else ; geo = {'Node' 'Bar' 'Tetrahedron' 'Hexahedron' 'Prism' 'Pyramid'} ; end
                    types = pkg.mesh.elements.AbstractElement.empty ;
                    for tt = 1:numel(geo) ; types(end+1) = eval(['pkg.mesh.elements.base.' geo{tt}]) ; end
                otherwise
                    error('Unsupported option for element types.') ;
            end
        end
    % Count nodes in each existing element
        nNodesInElems = this.Elems.nNodes ;
    % Check that types match node numbers
        valid = arrayfun(@(t)ismember(t.nNodes,nNodesInElems),types) ;
    % Assign
        this.Elems.Types = types(valid) ;
    end
end


%% MESH FACES & EDGES
properties (SetAccess = private)
    % Table of faces 
    Faces pkg.mesh.elements.ElementTable
    % Table of edges
    Edges pkg.mesh.elements.ElementTable
end
properties (Hidden)
    % Elem to features
    ElemFaces
    ElemEdges
end
methods % Compute the faces & edges only when needed
    function edges = get.Edges(this)
        if ~isempty(this.Elems) && isempty(this.Edges)
            [this.Edges,this.ElemEdges] = this.Elems.getTableOfUnique('Edges') ; 
        end
        edges = this.Edges ;
    end
    function faces = get.Faces(this)
        if ~isempty(this.Elems) && isempty(this.Faces) 
            [this.Faces,this.ElemFaces] = this.Elems.getTableOfUnique('Faces') ; 
        end
        faces = this.Faces ;
    end
end


%% COUNT FUNCTIONS
methods
    % Number of ...
    function N = nNodes(this) ; N = double(max(this.Elems.NodeIdx(:))) ; end
    function N = nElems(this) ; N = this.Elems.nElems ; end
    function N = nEdges(this) ; N = this.Edges.nElems ; end
    function N = nFaces(this) ; N = this.Faces.nElems ; end
    % Other ...
    function N = nMaxNodesByElem(this) ; N = this.Elems.nMaxNodesByElem ; end
    function N = nNodesInElem(this) ; N = this.Elems.nNodes ; end
end


%% MESH SPACE DIMENSION
properties (Dependent) ; nCoord ; end
methods
    % Allow to easily change the mesh space dimension (any integer>0)
    function set.nCoord(this,nCoord)
        this.X = [this.X.Values(:,1:min(this.nCoord,nCoord)) zeros(this.nNodes,nCoord-this.nCoord)] ;
    end
    function nCoord = get.nCoord(this) ; nCoord = size(this.X.Values,2) ; end
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
    % Elements to Edges [nEdges nElems]: edg(:) = M*elmt(:)
        M = this.Elems.contains(this.Faces) ;
    end
    function M = elem2edge(this,varargin)
    % Elements to Edges [nEdges nElems]: edg(:) = M*elmt(:)
        M = this.Elems.contains(this.Edges) ;
    end
    function M = face2edge(this,varargin)
    % Faces to Edges [nFaces nElems]: edg(:) = M*face(:)
        M = this.Faces.contains(this.Edges) ;
    end
end

%% BOUNDARIES
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


%% OPERATIONS ON CONNECTIVITIES (node/edge/face splitting, merging, etc)
methods
    function mesh = splitNodes(this,nod)
    % Split the mesh at given nodes
    % Create a new node for each attached element
        if islogical(nod) ; nod = find(nod) ; end
        nod = unique(nod(:)) ;
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        for nn = nod(:)'
            ind = ismember(mesh.Elems.NodeIdx,nn) ;
            newNodes = [nn mesh.nNodes + (1:sum(ind(:))-1)] ;
            mesh.Elems.NodeIdx(ind(:)) = newNodes(:) ;
            mesh.X.Values(newNodes,:) = repmat(mesh.X.Values(nn,:),[numel(newNodes) 1]) ;
        end
    end
    
    function mesh = splitEdges(this,edg)
    % Split the mesh at given edge indices
        if islogical(edg) ; edg = find(edg) ; end
        edg = unique(edg(:)) ;
    % Corresponding end nodes
        nodeIdx = this.Edges.NodeIdx(edg,:) ;
        nodeIdx = [nodeIdx(:,1) nodeIdx(sub2ind(size(nodeIdx),(1:size(nodeIdx,1))',sum(nodeIdx>0,2)))] ;
        nodes = unique(nodeIdx(:)) ;
    % Connectivities
        edg2nod = this.edge2node ;
        elm2edg = this.elem2edge ;
    % Browse over the nodes...
        elems = this.Elems ;
        x = this.X.Values ;
        for nn = 1:numel(nodes)
        % Edges linked to the node
            nodEdg = find(edg2nod(nodes(nn),:)) ;
        % End edges: do not share elements with the edges to split
            endEdg = sum(elm2edg(nodEdg,:)*elm2edg(edg,:)',2)==0 ;
        % Remove end edges
            normalEdg = nodEdg(~endEdg) ;
            normalEdg = setdiff(normalEdg,edg) ;
        % Create a new node for each normal edge
            newNodes = [nodes(nn) size(x,1)+(1:numel(normalEdg)-1)] ;
            x(newNodes,:) = repmat(x(nodes(nn),:),[numel(normalEdg) 1]) ;
        % Assign the new nodes to the normal edge's neightboring elements
            for ee = 1:numel(normalEdg)
                elmt = elm2edg(normalEdg(ee),:) ;
                idx = ismember(elems.NodeIdx,nodes(nn)) & elmt(:) ;
                elems.NodeIdx(idx) = newNodes(ee) ;
            end
        end
    % Assign the mesh
        if nargout==0 ; this.X = x ; this.Elems = elems ; 
        else ; mesh = pkg.mesh.Mesh('X',x,'Elems',elems) ; 
        end
    end
end


%% GEOMETRY
methods

    function C = circumCenters(this,table)
    % Return the circumcenters associated to the given list
    % By default, list is this.Elems
        if nargin<2 ; table = this.Elems ; end
        C = table.meanDataAtIndices(this.X.Values) ;
        C = reshape(C,table.nElems,this.nCoord) ;
    end
    
    function [ispl,normal,origin,frame] = isPlanar(this,tol)
    % Check for the mesh planarity
    % Also return the plane normal, origin 
    % and full frame [ivect icoord]
        origin = mean(this.X.Values,1) ;
        % If the mesh is trivially plane
            if this.nCoord<3 % 1D or 2D coordinates
                ispl = true ; 
                normal = [0 0 1] ;
                frame = eye(3) ;
                return ;
            end
        % Otherwise, try to find a common plane by finding the
        % principal coordinates variances
            if nargin<2 ; tol = 1e-9 ; end
            % Compute the coordinates covariance
                P = this.X.Values-origin ;
                [V,s] = eig(P'*P,'vector') ;
            % V is the frame basis, s the variances
                [s,ind] = sort(s,'descend') ;
                V = V(:,ind) ;
            % Test the planarity
                ispl = s(end)/s(1)<tol ;
            % Return objects
                frame = V' ;
                normal = frame(end,:) ;
    end
    
    function [frames,normals] = faceFrames(this,E)
    % Return the face frames of the mesh queried at local coordinates E
    % E: [nFaces nFaceDims==2]
    % frames: [nFaces 3 3] cat(3,x1,x2,normals)
    % normals: [nFaces 3]
        frames = NaN(this.Faces.nElems,3,3) ;
        if this.Faces.nElems==0 ; return ; end
        if nargin<2
            E = [0 0 ; this.Faces.Types.circumcenter] ;
            E = E(this.Faces.TypeIdx+1,:) ;
        end
        Xe = this.Faces.dataAtIndices(this.X.Values) ;
        for tt = 1:this.Faces.nTypes
            elemType = this.Faces.Types(tt) ;
            elemIdx = this.Faces.TypeIdx==tt ;
            J = elemType.evalJacobianAt(E(elemIdx,:)) ;
        end
        if nargin>1 ; normals = frames(:,:,3) ; end
    end
    
    function [frames,normals] = edgeFrames(this,E)
    % Return the edge frames of the mesh queried at local coordinates E
    % E: [nEdges nEdgeDims==1]
    % frames: [nEdges nCoord nCoord] cat(3,tangents,normals,(thirdvect))
    % normals: [nEdges nCoord]
        if this.Edges.nElems==0 ; return ; end
        if nargin<2
            E = [0 ; this.Edges.Types.circumcenter] ;
            E = E(this.Edges.TypeIdx+1,:) ;
        end
    % Compute the tangents
        tangents = NaN(this.Edges.nElems,this.nCoord) ;
        Xe = this.Edges.dataAtIndices(this.X.Values) ; % [nEdges nMaxNodesInEdges nEdgeDims==1]
        for tt = 1:this.Edges.nTypes
            elemType = this.Edges.Types(tt) ;
            elemIdx = this.Edges.TypeIdx==tt ;
            xe = Xe(elemIdx,1:elemType.nNodes,:) ; % [nE nNodesInElem nCoord]
            dN_de = elemType.evalJacobianAt(E(elemIdx,:)) ; % [nE nNodesInElem nEdgeDims==1]
            tangents(elemIdx,:) = permute(sum(permute(dN_de,[1 2 4 3]).*xe,2),[1 3 4 2]) ; % [nE nCoord nEdgeDims==1]
        end
    % Normalize
        tangents = tangents./sqrt(sum(tangents.^2,2)) ;
    % Rotate 90° about z axis for the normals
        if this.nCoord==3 ; normals = tangents*[0 -1 0 ; 1 0 0 ; 0 0 1] ;
        else ; normals = tangents*[0 -1 ; 1 0] ;
        end
    % Full Frames
        frames = cat(3,tangents,normals) ;
    % Cross product for the third vector
        if this.nCoord==3
            frames(:,:,3) = tangents(:,[2 3 1]).*normals(:,[3 1 2]) - tangents(:,[3 1 2]).*normals(:,[2 3 1]) ;
        end
    end

    function mesh = sortElems(this,dir)
    % Sort the element nodes of a planar mesh in clockwise 
    % or counter-clockwise(==trigo) direction
    % WARNING: only works with 1st order TRIS or QUADS !!!
        if ~this.isPlanar ; error('Cannot sort nodes of a non-planar mesh') ; end
        if nargin<2 ; dir = 'trigo' ; end
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        % Barycentered coordinates
            x = mesh.Elems.dataAtIndices(mesh.X.Values) ;
            x = x-mean(x,2,'omitnan') ;
        % Angle
            aaa = angle(x(:,:,1)+1i*x(:,:,2)) ;
        % Sort indices
            [~,indSort] = sort(aaa,2) ;
            switch dir
                case 'trigo'
                case 'counter-clockwise'
                case 'clockwise'
                    indSort = flip(indSort,2) ;
            end
        % To linear indices
            indSort = sub2ind(size(mesh.Elems.NodeIdx),repmat((1:mesh.nElems)',[1 mesh.nMaxNodesByElem]),indSort) ;
        % New element list
            mesh.Elems.NodeIdx = mesh.Elems.NodeIdx(indSort) ;
        % Re-set 1st quads node order
            Q4 = mesh.Elems.isP1Quad ; if any(Q4) ; mesh.Elems.NodeIdx(Q4,:) = mesh.Elems.NodeIdx(Q4,[1 2 4 3]) ; end
    end
    
    function sz = elemSize(this)
    % Compute the typical size of each element
    % 1D: length , 2D: area , 3D: volume
    % <TODO> Generalize this function
        sz = zeros(this.nElems,1) ;
        Xe = this.Elems.dataAtIndices(this.X.Values) ;
        isQ4 = this.Elems.isP1Quad ; if any(isQ4) ; Xe(isQ4,:,:) = Xe(isQ4,[1 2 4 3],:) ; end % Quad elements
        sz = polyarea(Xe(:,:,1),Xe(:,:,2),2) ;
%         for tt = 1:this.Elems.nTypes
%             elmtIdx = this.Elems.TypeIdx==tt ;
%             elmtType = this.Elems.Types(tt) ;
%             dN_de = elmtType.evalJacobianAt(elmtType.GaussIntegrationPoints) ; % [nIntPts nNodesInElmt nElmtDims]
%             dX_de = Xe(elmtIdx,1:elmtType.nNodes,:).*permute(dN_de,[4 2 5 3 1]) ; % [nElmt nNodesInElmt nCoord nElmtDims nIntPts]
%             dX_de = sum(dX_de,2) ; % [nElmt 1 nCoord nElmtDims nIntPts]
%             dX_de = sum(dX_de,4) ; % [nElmt 1 nCoord 1 nIntPts]
%             dX_de = prod(dX_de,3) ; % [nElmt 1 1 1 nIntPts]
%             sz(elmtIdx) = sum(permute(dX_de,[1 5 2 3 4]).*elmtType.GaussIntegrationWeights',2) ; % [nElmt 1]
%         end
    end
    
end


%% POINT LOCALIZATION & INTERPOLATION
methods
    
    function [M,pp,ee] = isInside(this,P,features,X,tol,bboxOnly)
    % Detect points P inside mesh (features=[]) or mesh features
    % Return a sparse logical matrix M(ii,jj)=true if P(ii) is (on/in)side feature(jj)
    % P = [nPts nCoord] double
    % features = [] (default) or ElementTable like Elems, Faces or Edges
    % X can be used to give other node coordinates
    % tol is the tolerance
    % bboxOnly can be used to check if a point is inside an element bbox
    % (skip the localization+insideelmt steps)
        if nargin<3 ; features = [] ; end
        if nargin<4 ; X = this.X.Values ; end
        if nargin<5 ; tol = this.defaultTolerance(X) ; end
        if nargin<6 ; bboxOnly = false ; end
        nPts = size(P,1) ;
    % Points in the global mesh bounding box ?
        globalBBox = [min(X,[],1) ; max(X,[],1)] ;
        in = all( P>=globalBBox(1,:)-tol & P<=globalBBox(2,:)+tol ,2) ;
        pp = find(in) ; ee = [] ;
    % If no features are given, return the current values
        if isempty(features) ; M = sparse(pp,pp*0+1,true,nPts,1) ; return ; end
    % If no points are inside the mesh Bounding box, return an empty matrix
        if isempty(pp) ; M = sparse(nPts,features.nElems) ; return ; end
    % Find the closest nodes 
        [nn,dist] = this.closestNode(P(pp,:),X) ;
        e2n = sparse(features) ;
        [ppp,ee] = find(e2n(nn,:)) ;
        pp = pp(ppp) ;
    % Points in element bounding box ?
        [pp,ee] = inElmtBBox(this,P,features,X,tol,pp,ee) ;
        if bboxOnly ; M = logical(sparse(pp,ee,1,nPts,features.nElems)) ; return ; end
    % Build the matrix
        M = logical(sparse(pp,ee,1,nPts,features.nElems)) ;
    end
    
    function [nn,dist] = closestNode(this,P,X)
    % Return the closest node indices
    % P = [nPts nCoord] or [nPts 1 nCoord]
    % P = [nNodes nCoord] or [nNodes 1 nCoord]
    % nn = [nPts 1]
        if nargin<3 ; X = this.X.Values ; end
        [nn,sqdist] = pkg.data.closestPoint(P,X) ;
        if nargout>1 ; dist = sqrt(sqdist) ; end
    end
    
    function [pp,ee] = inElmtBBox(this,P,features,X,tol,ip,ie)
    % Return indices pp and ee corresponding to P(pp) inside elmt(ee)
    % P = [nPts nCoord]
    % X = [nNodes nCoord]
    % By default, each point is tested against each element of the mesh
    % Optionally, (ip ie) [nTest 1] can be provided to test specific point with
    % specific elements (see below)
        if nargin<3 ; features = this.Elems ; end
        if nargin<4 ; X = this.X.Values ; end
        if nargin<5 ; tol = this.defaultTolerance(X) ; end
        nPts = size(P,1) ;
    % Put the coordinates to the third dimension
        P = reshape(P,nPts,1,[]) ; % [nPts 1 nCoord]
        if nargin>5 ; P = P(ip,:,:) ; end % [nTest 1 nCoord]
    % Element node coordinates
        Xe = features.dataAtIndices(X) ; % [nElems nMaxNodesByElem nCoord]
        if nargin>5 % test specific points vs specific elements
            Xe = Xe(ie,:,:) ; % [nTest nMaxNodesByElem nCoord] 
            elmtBBox = [min(Xe,[],2) max(Xe,[],2)] ; % [nTest 2 nCoord]
            in = all(P>=elmtBBox(:,1,:)-tol & P<=elmtBBox(:,2,:)+tol,3) ; % [nTest 1]
            ii = find(in) ;
            pp = ip(ii) ; ee = ie(ii) ;
        else % test all points vs all elements
            elmtBBox = [min(Xe,[],2) max(Xe,[],2)] ; % [nElems 2 nCoord]
            elmtBBox = permute(elmtBBox,[2 1 3]) ; % [2 nElems nCoord]
            in = all(P>=elmtBBox(1,:,:)-tol & P<=elmtBBox(2,:,:)+tol,3) ; % [nPts nElems]
            [pp,ee] = find(in) ;
        end
    end
    
    function [E,ie] = localize(this,P,features,extrap,X,tol,ip,ie)
    % Localize points P on the mesh. 
    % Returns the local coordinates corresponding to the feature's closest point
    % Find E = argmin(norm(features.evalAt(E)*X-P))
    % P = [nPts nCoord]
    % features = ElementTable: mesh Elems (default), Faces or Edges
    % extrap = return point localization even outside the features (else E(outside,:) = NaN)
    % extrap is by default false if (ip,ie) are not provided
    % X = [nNodes nCoord] (mesh.X by default)
    % tol: convergence tolerance
    % ip,ie = [nLocal 1] optionnal: restricts the test between P(ip) and features(ie)
    % E = [nPts nElems nMaxElemDim] or [nLocal nMaxElemDim]
        if nargin<3 ; features = this.Elems ; end
        if nargin<4 ; extrap = nargin<7 ; end
        if nargin<5 ; X = this.X.Values ; end
        if nargin<6 ; tol = this.defaultTolerance(X) ; end
        if nargin==7 ; error('Point AND element indices couples (ip,ie) must be provided') ; end
        itMax = 10 ;
        debug = false ;
        if debug ; pl0 = plot(P(:,1),P(:,2),'or') ; pl = plot(P(:,1)*NaN,P(:,2)*NaN,'.b') ; end
    % Localization mode
        if nargin>7 % Do what is asked, all inputs are available
            mode = 'custom' ;
        elseif extrap % Localize each point in each mesh feature
            mode = 'allInAll' ;
        else % Localize each point in the closest mesh features (bbox only)
            mode = 'inBBox' ;
        end
    % If extrapolation is allowed, localize each point in each mesh feature
        nPts = size(P,1) ; nElems = features.nElems ;
        switch mode
            case 'custom' % Do what is asked, all inputs are available
            case 'allInAll' % Do what is asked, all inputs are available
                ip = repmat(1:nPts,[1 nElems]) ; 
                ie = repelem(1:nElems,nPts*ones(1,nElems)) ; 
            case 'inBBox'
                [~,ip,ie] = isInside(this,P,features,X,tol,true) ; % true for bboxOnly
        end
        nLocal = numel(ip) ;
    % Reshape the list of points
        P = P(ip,:)' ; % [nCoord nLocal]
    % Coordinates of element nodes
        Xe = features.dataAtIndices(X) ; % [features.nElems features.nMaxNodesByElem nCoord]
    % Initialize the output
        E = NaN(nLocal,max([features.Types.nDims])) ;
    % For each element type
        for typeIdx = 1:features.nTypes
            elmtType = features.Types(typeIdx) ; % The type of feature to deal with
            % Find the corresponding element indices
                elmtIdx = find(features.TypeIdx==typeIdx) ; % The indices in the feature list
                [isType,ii] = ismember(ie,elmtIdx) ; % Keep only elements of interest
                elmtIdx = elmtIdx(ii(isType)) ; % The element indices we will have to deal with
                nE = length(elmtIdx) ;
            % Retrieve the element node coordinates
                xe = Xe(elmtIdx,1:elmtType.nNodes,:) ; % [nE nElmtNodes nCoord]
                xe = permute(xe,[3 1 2]) ; % [nCoord nE nElmtNodes]
            % Intialize the algorithm
                e0 = mean(elmtType.NodeLocalCoordinates,1) ;
                if norm(reshape(elmtType.evalHessianAt(e0),[],1))<tol | 1 ; methodOrder = 1 ;
                else ; methodOrder = 2 ; end
                e = repmat(e0,[nE 1]) ; % [nE nElmtNodes]
            % Loop until norm(x-P)<tol
                notConverged = NaN ; it = 0 ;
                while ~isempty(notConverged) && it<itMax
                    it = it + 1 ;
                % Shape function evalutation (+derivatives)
                    N = elmtType.evalAt(e) ; % [nE nElmtNodes]
                    dN_de = elmtType.evalJacobianAt(e) ; % [nE nElmtNodes nElmtDims]
                    if methodOrder>1 ; d2N_de2 = elmtType.evalHessianAt(e) ; end % [nE nElmtNodes nElmtDims nElmtDims]
                % Function guess evaluation (+derivatives)
                    x = sum(permute(N,[3 1 2]).*xe,3) ; % [nCoord nE]
                    if debug ; pl.XData = x(1,:) ; pl.YData = x(2,:) ; drawnow ; end
                    dx_de = sum(permute(dN_de,[4 1 2 3]).*xe,3) ; % [nCoord nE 1 nElmtDims]
                    if methodOrder>1 ; d2x_de2 = sum(permute(d2N_de2,[5 1 2 3 4]).*xe,3) ; end % [nCoord nE 1 nElmtDims nElmtDims]
                    dx_de = permute(dx_de,[1 4 2 3]) ; % [nCoord nElmtDims nE]
                    if methodOrder>1 ; d2x_de2 = permute(d2x_de2,[1 4 5 2 3]) ; end % [nCoord nElmtDims nElmtDims nE]
                % Residue check
                    res = x-P(:,isType) ; % [nCoord nE]
                    notConverged = find(sum(res.^2,1)>tol^2) ;
                    if isempty(notConverged) ; break ; end
                % Update points where it needs to be
                    for iii = notConverged
                        if methodOrder==1
                            de =  - dx_de(:,:,iii)\res(:,iii) ;
                        else
                            J = dx_de(:,:,iii) ;
                            H = J'*J + squeeze(sum(res(:,iii).*d2x_de2(:,:,:,iii),1)) ;
                            de =  - H\(J'*res(:,iii)) ;
                        end
                        e(iii,:) = e(iii,:) + de.' ;
                    end
                end
            % Assign
                E(isType,:) = e ;
            % If not extrapolating, check that E correspond to local
            % coordinate INSIDE the element
                if ~extrap
                    outside = ~elmtType.isInside(E(isType,:)) ;
                    E(isType & outside,:) = NaN ;
                end
        end
    % Reshape if needed
        if strcmp(mode,'allInAll') ; E = reshape(E,nPts,nElems,[]) ; end
        if debug ; delete(pl0) ; delete(pl) ; end
    end
    
end

%% MESH OPERATIONS (MOTION, BOOLEAN)
methods 

    function mesh = plus(this,mesh2)
    % Addition of two meshes (NEED CLEANUP TO WELD THE MESHES!)
        mesh = copy(this) ; % Always take a copy
        % New mesh coordinates values
            nNodes = this.nNodes + mesh2.nNodes ;
            x = NaN(nNodes,max(this.nCoord,mesh2.nCoord)) ;
            x(1:this.nNodes,1:this.nCoord) = this.X.Values ;
            x(this.nNodes+1:end,1:mesh2.nCoord) = mesh2.X.Values ;
        % New element table
            elems = [this.Elems mesh2.Elems] ;
            elems.NodeIdx(this.nElems+1:end,:) = elems.NodeIdx(this.nElems+1:end,:) + this.nNodes ;
        % Return/Modify the mesh
            mesh.X = x ; 
            mesh.Elems = elems ;
    end
    
    function meshes = splitParts(this)
    % Split the independent parts of a mesh into several-meshes
        e2n = this.elem2node ;
        e2e = e2n'*e2n ;
        meshes = pkg.mesh.Mesh.empty ;
        elems = this.Elems ;
        while ~isempty(e2e)
            belong = e2e\sparse(1,1,1,size(e2e,1),1)~=0 ;
            belongElems = pkg.mesh.elements.ElementTable('Types',this.Elems.Types,'Indices',this.Elems.Indices(belong,:)) ;
            meshes(end+1) = pkg.mesh.Mesh('X',this.X.Values,'Elems',belongElems) ;
            elems.Indices = elems.Indices(~belong,:) ;
            e2e = e2e(~belong,~belong) ;
        end
    end
    
    function mesh = merge(this,mesh2,tol)
    % Merge two meshes (WIP) [need mesh cleanup after that]
        if nargin<3 ; tol = this.defaultTolerance ; end
    % Find the closest nodes
        % Reduce to bounding boxes
            bboxThis = [min(this.X.Values,[],1) ; max(this.X.Values,[],1)] + 1.2*[-1 ; 1]*tol ;
            bboxMesh2 = [min(mesh2.X.Values,[],1) ; max(mesh2.X.Values,[],1)] + 1.2*[-1 ; 1]*tol ;
            thisNodesInMesh2 = find(all(this.X.Values>=bboxMesh2(1,:) & this.X.Values<=bboxMesh2(2,:),2)) ;
            mesh2NodesInThis = find(all(mesh2.X.Values>=bboxThis(1,:) & mesh2.X.Values<=bboxThis(2,:),2)) ;
        % Compute the distances
            sqDist = sum((permute(mesh2.X.Values(mesh2NodesInThis,:),[1 3 2])-permute(this.X.Values(thisNodesInMesh2,:),[3 1 2])).^2,3) ;
            [nodMesh2,nodThis] = find(sqDist<=tol.^2) ; % all pairs of close nodes
            sqDist = sqDist(sub2ind(size(sqDist),nodMesh2,nodThis)) ;
        % Convert to global numerotation (including out-of-bbox)
            nodThis = thisNodesInMesh2(nodThis) ;
            nodMesh2 = mesh2NodesInThis(nodMesh2) ;
    % Do not merge multiple nodes with one node (!!!)
        % Sort by increasing distance
            [~,indSort] = sort(sqDist,'ascend') ;
            nodMesh2 = nodMesh2(indSort) ; 
            nodThis = nodThis(indSort) ;
        % Take only the closest node that is not already taken
            [nodThis,nn] = unique(nodThis,'stable') ;
            nodMesh2 = nodMesh2(nn) ;
            [nodMesh2,nn] = unique(nodMesh2,'stable') ;
            nodThis = nodThis(nn) ;
    % Merge Nodes
        newNodes = mean(cat(3,this.X.Values(nodThis,:),mesh2.X.Values(nodMesh2,:)),3) ;
        % Complete list
            x = [this.X.Values ; mesh2.X.Values] ;
            newNodeIdx = 1:size(x,1) ;
        % Replace merged nodes
            x(nodThis,:) = newNodes ;
            x(this.nNodes + nodMesh2,:) = newNodes ;
        % Cull duplicates
            newNodeIdx(this.nNodes + nodMesh2) = nodThis ;
            newNodeIdx(setdiff(1:numel(newNodeIdx),this.nNodes+nodMesh2)) = 1:numel(newNodeIdx)-numel(nodMesh2) ;
            x(this.nNodes + nodMesh2,:) = [] ;
    % Merge elements
        elems = [this.Elems ; mesh2.Elems] ; 
        elems.NodeIdx(this.nElems+1:end,:) = elems.NodeIdx(this.nElems+1:end,:) + this.nNodes ;
        newNodeIdx = [0 newNodeIdx] ;
        elems.NodeIdx = newNodeIdx(elems.NodeIdx+1) ;
    % Return the new mesh
        mesh = pkg.mesh.Mesh('X',x,'Elems',elems) ;
    end

    function mesh = applyTransform(this,F,v)
    % Apply a geometrical transformation to the mesh
    % F is the transfo: [nCoord nCoord N] array,
    % v is the translation vector; has to be of size [N nCoord]
    % with N>=1 to replicate the mesh on different transfos.
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        if nargin<3 ; v = zeros(1,size(F,2)) ; end
        % Force formatting
            F = F(:,:,:) ; v = permute(v(:,:),[3 2 1]) ; % [1 nCoord N]
            F = F + v*0 ; v = v + F(1,:,:)*0 ;
            [nC,~,N] = size(F) ;
        % Fix incompatible dimensions
            dC = mesh.nCoord-nC ;
            if dC>0
                O = zeros(nC,dC,N) ;
                F = [F O ; permute(O,[2 1 3]) repmat(eye(dC),[1 1 N])] ;
                v = [v O(1,:,:)] ;
            elseif dC<0
                mesh.nCoord = nC ;
            end
        % New element list if needed
            elems = repmat(mesh.Elems,[N 1]) ;
            elems.NodeIdx = elems.NodeIdx + uint32(mesh.nNodes*kron((0:N-1)',ones(mesh.nElems,1))) ;
        % Add a fourth dimension to stack replicated meshes
            F = permute(F,[4 1 2 3]) ; % [1 nCoord nCoord N]
            v = permute(v,[4 1 2 3]) ; % [1 1 nCoord N]
        % Apply Transform
            x = sum(mesh.X.Values.*F,2) + v ; % [nNodes 1 nCoord N]
        % Reshape
            x = permute(x,[1 4 3 2]) ; % [nNodes N nCoord]
            x = reshape(x,[this.nNodes*N this.nCoord]) ; % [nNodes*N nCoord]
        % Set
            mesh.X = x ;
            mesh.Elems = elems ;
    end

    function mesh = move(this,T)
    % Move the mesh by a translation vector T
    % T has to be of size [N nCoord] with N>=1 to replicate the mesh
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        % Format the translation vector
            T = squeeze(T) ;
        % Apply transform
            mesh.applyTransform(eye(size(T,2)),T) ;
    end

    function mesh = rotate(this,THETA,CEN,AX)
    % Rotate the mesh by
    % the angle THETA: [N 1]  
    % around the point P: [1 nCoord] (default [0 0 0])
    % with respect to the axis AX: [1 3] (default [0 0 1], not implemented)
    % with N>=1 to replicate the mesh
        if nargin<3 ; CEN = zeros(1,mesh.nCoord) ; end
        if nargin<4 ; AX = [0 0 1] ; end
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        % Move the mesh
            mesh.move(-CEN) ;
        % Compute the transformation matrix
            THETA = reshape(THETA(:),1,1,[]) ;
            F = [cos(THETA) sin(THETA); -sin(THETA) cos(THETA)] ;
            if mesh.nCoord>2 ; F(end+1,end+1,:) = 1 ; end
        % Rotate
            mesh.applyTransform(F) ;
        % Re-move the mesh
            mesh.move(CEN) ;
    end
    
    function mesh = replicate(this,x,fill,loop)
    % Replicate a mesh with the given 3D-array coordinates
    % X can either be: 
    %   - array [this.nNodes nCoord nRep] of new coordinates
    %   - int. scalar: number of replicates
    % fill (bool) if the replicates are joined to a sweeped mesh
    % loop (bool) close a filled mesh (join the first replicate with the
    % last)
        if nargin<2 ; error('Not enough arguments') ; end
        if numel(x)==1 ; x = repmat(this.X.Values,[1 1 x]) ; end
        if size(x,1)~=size(this.X.Values,1) ; error('Wrong format for input X.') ; end
        nRep = size(x,3) ;
        if nargin<3 ; fill = nRep>1 ; end
        if nargin<4 ; loop = false ; end
    % New node coordinates
        x = reshape(permute(x,[1 3 2]),this.nNodes*nRep,[]) ;
    % New element indices
        indices = double(this.Elems.Indices) ; % convert for following operations
        indices = repmat(indices,[1 1 nRep]) ; % replicate elements
        indices(:,2:end,:) = indices(:,2:end,:) + reshape(0:nRep-1,[1 1 nRep])*this.nNodes ; % update node indices (not TypeIdx)
        indices = indices.*logical(indices(:,:,1)) ; % re-set invalid indices to 0
    % Fill ?
        types = this.Elems.Types ;
        if fill % need to return extruded elements
            if loop ; indRep = [1:nRep ; [2:nRep,1]] ;
            else ; indRep = [1:nRep-1 ; 2:nRep] ;
            end
            nElems = this.nElems*size(indRep,2) ; % final number of elements
            nMaxNodes = this.nMaxNodesByElem ; % current maximum number of nodes
            indices = permute(indices,[1 3 2]) ; % [this.nElems nRep nMaxNodes+1]
            indices = [ reshape(indices(:,1:end-1,1),nElems,1) ...
                        reshape(indices(:,indRep(1,:),2:end),nElems,nMaxNodes) ...
                        reshape(indices(:,indRep(2,:),2:end),nElems,nMaxNodes) ...
                        ] ;
            for tt = 1:numel(types)
            % Valid nodes (for heterogeeous element types...
                validNodes = (1:types(tt).nNodes)'+[0 1]*nMaxNodes ;
            % Depending on the initial element type...
                switch class(types(tt))
                    case 'pkg.mesh.elements.base.Bar'
                        types(tt) = pkg.mesh.elements.base.Quadrangle ;
                        validNodes(:,2) = flip(validNodes(:,2)) ; % quad numerotation...
                    case 'pkg.mesh.elements.base.Triangle'
                        types(tt) = pkg.mesh.elements.base.Prism ;
                    case 'pkg.mesh.elements.base.Quadrangle'
                        types(tt) = pkg.mesh.elements.base.Hexahedron ;
                    otherwise
                        error(['Cannot extrude a ' class(types(tt)) ' !']) ;
                end
            % Change indices
                elmtIdx = indices(:,1)==tt ;
                indices(elmtIdx,:) = [  indices(elmtIdx,1) ...
                                        indices(elmtIdx,validNodes(:)+1) ...
                                        zeros(sum(elmtIdx),size(indices,2)-types(tt).nNodes-1)...
                                        ] ;
            end
        else
            indices = indices(:,:) ; % to 2d array of indices
        end
        elems = pkg.mesh.elements.ElementTable('Types',types,'Indices',indices) ;
    % Assign
        if nargout==0 ; this.X = x ; this.Elems = elems ; 
        else ; mesh = pkg.mesh.Mesh('X',x,'Elems',elems) ; 
        end
    end
    
    function mesh = offset(this,DIST,fill)
    % Offset 2D/3D mesh edges/faces
    % DIST is the offset distances (can be negative)
    % fill: (bool) set true when you want to create..
    % ..a surface/volume from edges/faces
    % fill: bool. false: just offset the mesh; true: fill with elements
        if any(this.Elems.Types.nDims>3) ;  error('Cannot offset 3D elements') ; end
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
    % Process inputs
        if numel(DIST)==0 ; return ; end % return the same mesh
        if nargin<3 ; fill = numel(DIST)>1 ; end % fill by default
        if fill && numel(DIST)==1 ; DIST = [0 ; DIST(:)] ; end % force the volume to be valid
        DIST = sort(DIST(:),'ascend') ; 
        nOfst = numel(DIST) ;
    % Operate
        % PLANE WIRE-ONLY MESH
        if this.isPlanar && all(this.Elems.Types.nDims==1)
            % Compute node normals
                [~,edgeNormals] = edgeFrames(this) ;
                nodeNormals = this.edge2node('mean')*edgeNormals ;
            % Offset Nodes
                offset = nodeNormals.*reshape(DIST,[1 1 nOfst]) ;
                x = this.X.Values + offset ;
        % SURFACE-ONLY MESH
        elseif all(this.Elems.Types.nDims==2) 
            error('Surface offset is not implemented yet') ;
        end
    % Replicate elements with the given coordinates array
        mesh.replicate(x,fill) ;
    end
    
    function mesh = sweep(this,CRV,mode)
    % Create a sweep mesh geometry from a mesh and a curve
    % CRV: [nPts nCoord]
    % mode: 'angle' or 'parallel' (default)
    %   - 'angle': the replicated meshes remain at the same angle w.r.t the curve
    %   - 'parallel': all nodes are swept by the same curve (parallel edges)
        if any(this.Elems.Types.nDims>3) ; error('Cannot offset 3D elements') ; end
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
    % Process inputs
        if nargin<3 ; mode = 'parallel' ; end
    % Initialize
        x = this.X.Values ;
    % Match Curve & nodes dimensions
        nC = max(size(x,2),size(CRV,2)) ; % target nCoord
        x = [x zeros(size(x,1),nC-size(x,2))] ;
        CRV = [CRV zeros(size(CRV,1),nC-size(CRV,2))] ; % [nPts nC]
    % Sweep nodes
        u0 = permute(CRV-CRV(1,:),[3 2 1]) ; % [1 nC nPts] ;
        switch mode
            case 'parallel' % x = x + u0
                x = x + u0 ; % [nNodes nC nPts]
            case 'angle' % x = x + u0 + vectprod(R,z) 
                % Tangent vector [nPts nC]
                    tang = diff(CRV,1,1) ; % segments
                    tang = [3*tang(1,:)-tang(2,:) ; tang(1:end-1,:)+tang(2:end,:) ; 3*tang(end,:)-tang(end-1,:)] ; % node tangents
                    tang = tang./sqrt(sum(tang.^2,2)) ; % normalize
                % Rotation vector R = vectprod(t_{n},t_{n+1})
                    R = tang(1,[2 3 1]).*tang(:,[3 1 2]) - tang(1,[3 1 2]).*tang(:,[2 3 1]) ; % [nPts nC]
                    R = permute(R,[3 2 1]) ; % [1 nC nPts]
                % Apply
                    z = x-CRV(1,:) ; % relative coordinates
                    r = [ R(:,2,:).*z(:,3,:) - R(:,3,:).*z(:,2,:) ...
                          R(:,3,:).*z(:,1,:) - R(:,1,:).*z(:,3,:) ...
                          R(:,1,:).*z(:,2,:) - R(:,2,:).*z(:,1,:) ] ;
                    x = x + u0 + r ;
            otherwise
                error('Unknown sweeping mode')
        end
    % Replicate elements given the new coordinates array
        mesh.replicate(x,true) ;
    end
    
    function mesh = extrude(this,VEC,DIV)
    % Extrude the mesh along the vector VEC
    % (opt.) use DIV (int. scalar) to subdivide the mesh
    % VEC: [1 nCoord]
        if nargin<3 ; DIV = 1 ; end
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        VEC = VEC(:)' ;
        CRV = interp1(VEC.*[0;1],linspace(1,2,DIV+1)) ;
        x = mesh.X.Values + permute(CRV,[3 2 1]) ;
        mesh.replicate(x,true) ; 
    end
    
    function mesh = revolution(this,AX,PT,ANG,DIV,loop) 
    % Revolve a mesh around AX centered on PT with angle ranges ANG
    % AX: revolution axis [1 nCoord] (default [0 0 1])
    % PT: origin point [1 nCoord] (default [0 0 0])
    % ANG: revolution angle range (default [0 2/pi])
    %   - scalar A: ANG = [0 A]
    %   - range A: ANG = [Amin Amax]
    %   - vector [1 nA]: gives the subdivision
    % DIV: uniform subdivision
        if nargin<2 ; AX = [0 0 1] ; end 
        AX = AX(:)' ; AX = AX/norm(AX) ;
        if nargin<3 ; PT = [0 0 0] ; else ; PT = PT(:)' ; end
        if nargin<4 ; ANG = [0 2*pi] ; else ; ANG = ANG(:)' ; end
        if numel(ANG)==1 ; ANG = sort([0 ANG]) ; end
        if nargin<5 % We have to set DIV ?
            ANG = linspace(min(ANG(:)),max(ANG(:)),100) ; 
        else % DIV is given
            ANG = linspace(min(ANG(:)),max(ANG(:)),DIV+1) ;
        end
        DIV = numel(ANG)-1 ;
        if nargin<6 ; loop =  abs(ANG(end)-ANG(1)-2*pi)<100*eps ; end
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
    % Rotation center point P
        C = mean(mesh.X.Values,1,'omitnan') ; % current mesh barycenter
        P = PT + AX*sum((C-PT).*AX,2) ; % closest point on the revolution axis
    % Rotation matrices
        R = pkg.math.rotmat(-ANG(:),AX) ;
    % New coordinates
        x = pkg.math.innerprod(mesh.X.Values-P,R,2,1) + P ;
    % Apply
        mesh.replicate(x,true,loop) ;
    end   
    
end


%% MESH SLICING, CUTTING, SPLITTING,... BY A LEVEL SET
methods
    
    function sBool = lvlSetSign(this,fcnX,tol,X)
    % Evaluate a levelset f = fcn(X) on mesh nodes X and return a signed
    % boolean: -1 (inside), 0 (on), 1 (outside)
        if nargin<3 ; tol = this.defaultTolerance ; end
        if nargin<4 ; X = this.X.Values ; end
        if ~isnumeric(fcnX) ; fcnX = fcnX(X) ; end % fcnX can be a function handle..
        %sBool = sign(floor(abs(fcnX)/tol).*sign(fcnX)) ; % With a tolerance on the 'ON' case (fcn~0)
        sBool = sign(fcnX) ; % Without a tolerance
    end
    
    function [fBool,eBool,sBool] = featureLevelSetSign(this,fcnX,features,tol,X)
    % Same as above but at the level of features: edges/faces/elems
    % fBool: [nFeat 1] feature level
    % eBool: [nFeat nNaxNodesByFeature] feature node's level
    % sBool: [nNodes 1] node level
        if nargin<3 ; features = this.Elems ; end
        if nargin<4 ; tol = this.defaultTolerance ; end
        if nargin<5 ; X = this.X.Values ; end
        if ~isnumeric(fcnX) ; fcnX = fcnX(X) ; end % fcnX can be a function handle..
        sBool = lvlSetSign(this,fcnX,tol,X) ;
    % Spread in features
        eBool = features.dataAtIndices(sBool) ;
    % Set
        fBool = zeros(features.nElems,1) ; % on the level set by default (=0)
        fBool(all(eBool==-1 | isnan(eBool),2)) = -1 ; % all nodes inside
        fBool(all(eBool==1 | isnan(eBool),2)) = 1 ; % all nodes outside
    end
    
    function [Xc,t] = edgCrossLvlSet(this,fcnX,edg,tol,X)
    % Return the parameters t corresponding to the location where the edge
    % cross a levelset (fcn changes of sign)
        if nargin<3 ; edg = 1:this.nEdges ; end
        if nargin<4 ; tol = this.defaultTolerance ; end
        if nargin<5 ; X = this.X.Values ; end
        if isempty(edg) || ~any(edg) ; t = [] ; Xc = [] ; return ; end
        edges = this.Edges.subpart(edg) ;
        if ~isnumeric(fcnX) ; fcnX = fcnX(X) ; end
    % Get the intersection parameter <TODO> higher-order interp with tol..
        fcnEdg = edges.dataAtIndices(fcnX) ;
        t = fcnEdg(:,1)./(fcnEdg(:,1)-fcnEdg(:,2)) ; % (linear interpolation !)
    % Get intersection points
        Xe = edges.dataAtIndices(X) ;
        Xc = Xe(:,1,:).*(1-t) + Xe(:,2,:).*t ;
        Xc = permute(Xc,[1 3 2]) ;
    end
    
    function mesh = slice(this,fcn,tol,X)
    % Slice of the mesh at locations X where the sign of fcn(X) changes
        if nargin<3 ; tol = this.defaultTolerance ; end
        if nargin<4 ; X = this.X.Values ; end
    % Evaluate the levelset on nodes
        fcnX = fcn(X) ;
    % Get which elements need to be sliced
        [fBool,eBool,sBool] = featureLevelSetSign(this,fcnX,this.Elems,tol,X) ;
        indElemsToSlice = find(fBool==0) ;
        eBool = eBool(indElemsToSlice,:) ;
        elemsToSlice = this.Elems.subpart(indElemsToSlice) ;
    % Edges: we need to retrieve the connectivity...
        [edg,elmt,elmtIdx] = find(this.ElemEdges(:,indElemsToSlice)) ;
        globalEdgIdx = full(sparse(elmt,elmtIdx,edg)) ;
        elemEdges = pkg.mesh.elements.ElementTable('NodeIdx',globalEdgIdx) ;
    % Slice elements
        sliceElems = pkg.mesh.elements.ElementTable ;
        for tt = 1:elemsToSlice.nTypes
            % Retrieve concerned elements
                elmtType = elemsToSlice.Types(tt) ;
                elmtIdx = elemsToSlice.TypeIdx==tt ;
            % Slice the corresponding elements
                [newElems,elmtIdx] = elmtType.slice(eBool(elmtIdx,1:elmtType.nNodes)) ;
            % From local node indices to global node indices... (tidy part..)
            if ~isempty(elmtIdx)
                % Separate node & edge indices
                    edgIdx = newElems.NodeIdx.*uint32(newElems.NodeIdx>elmtType.nNodes) ;
                    nodIdx = newElems.NodeIdx-edgIdx ;
                % Go global
                    nodIdx = elemsToSlice.globalize(nodIdx,elmtIdx) ;
                    edgIdx = elemEdges.globalize(max(edgIdx-elmtType.nNodes,0),elmtIdx) ;
                    edgIdx(edgIdx>0) = edgIdx(edgIdx>0) + this.nNodes ;
                % Assign new elements
                    newElems.NodeIdx = nodIdx + edgIdx ;
                % Add to the list
                    sliceElems = [sliceElems newElems] ;
            end
        end
    % Clean the element list (may be duplicated elements)
        sliceElems = clean(sliceElems) ;
    % Unique (and shorter..) list of nodes/edges indices
        [elmtIdx,nodePos,idx] = find(sliceElems.NodeIdx) ;
        isNod = idx<=this.nNodes ; isEdg = ~isNod ; % imaginary indices denote edge indices
        [uNod,na,nc] = unique(idx(isNod)) ;
        [uEdg,ea,ec] = unique(idx(isEdg)-this.nNodes) ;
    % Keep nodes
        Xn = X(uNod,:) ;
        idx(isNod) = nc ;
    % Create new nodes at edge intersections
        Xc = edgCrossLvlSet(this,fcn,uEdg,tol,X) ;
        idx(isEdg) = ec + numel(uNod) ;
    % Finally..
        sliceElems.NodeIdx = full(sparse(elmtIdx,nodePos,double(idx))) ;
        Xs = [Xn ; Xc] ;
    % Slice mesh
        mesh = pkg.mesh.Mesh('X',Xs,'Elems',sliceElems) ;
    end
    
    function mesh = cut(this,fcn)
    % Return a cut of the mesh so that for any X in mesh, fcn(X)<=0
    end
    
    function mesh = split(this,fcn)
    % Split the mesh at locations where the sign of fcn(X) changes
    end
end


%% MESH CLEANUP
methods

    function tol = defaultTolerance(this,X)
    % Return the default tolerance used in cleanup operations
    % Optional X [nNodes nCoord] can be used for custom mesh coordinates
        if nargin<2 ; X = this.X.Values; end
        tol = norm(range(X,1))*1e-9 ;
    end

    function mesh = clean(this,tol)
    % Apply the following cleanup functions to the mesh
        if nargin<2 ; tol = this.defaultTolerance ; end
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        mesh.cullInvalid ;
        mesh.cullUnused ;
        mesh.cullDuplicates(tol) ;
        this.Elems = clean(this.Elems) ;
    end

    function mesh = removeNodes(this,remove)
    % Remove nodes in the list and clean the mesh
    % nod can be logical or indexes
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        if ~islogical(remove) % convert to logical
            remove = full(sparse(remove(:),ones(numel(remove),1),true,mesh.nNodes,1)) ; 
        end
        if ~any(remove) ; return ; end
        mesh.X = mesh.X.Values(~remove,:) ;
        mesh.Elems.NodeIdx = this.Elems.dataAtIndices(cumsum(~remove(:))) ;
    end

    function mesh = removeElems(this,remove)
    % Remove elems in the list and clean the mesh
    % nod can be logical or indexes
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        if ~islogical(remove) % convert to logical
            remove = full(sparse(remove(:),ones(numel(remove),1),true,mesh.nElems,1)) ; 
        end
        if ~any(remove) ; return ; end
        mesh.Elems.Indices = mesh.Elems.Indices(~remove,:) ;
        mesh.cullUnused ;
    end

    function iNodes = invalidNodes(this)
    % Detect invalid nodes (Inf or NaN)
        iNodes = any( isinf(this.X.Values) | isnan(this.X.Values) ,2) ;
    end

    function iElems = invalidElems(this)
    % Detect invalid elements (some nodes are invalid)
        iElems = any(this.Elems.NodeIdx>this.nNodes,2) ; % At least one node index too large
        iElems = iElems | all(this.Elems.NodeIdx<=0,2) ; % All invalid indices
        % Contains invalid nodes
        iNodes = invalidNodes(this) ;
        iElems = iElems | any(this.Elems.dataAtIndices(iNodes),2) ;
    end

    function mesh = cullInvalid(this)
    % Cull invalid nodes and elements
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        mesh.removeElems(mesh.invalidElems) ;   
        mesh.removeNodes(mesh.invalidNodes) ;         
    end

    function uNodes = unusedNodes(this)
    % Return true for nodes not part of element list
        uNodes = ~ismember(1:this.nNodes,this.Elems.NodeIdx(:)) ;
    end

    function mesh = cullUnused(this)
    % Cull unused nodes
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        mesh.removeNodes(mesh.unusedNodes) ;
    end

    function [mesh,nodesMoved] = cullDuplicates(this,tol)
    % Cull duplicate nodes in the mesh and corresponding elements
        if nargin<2 ; tol = this.defaultTolerance ; end
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
    % Unique node indices list
        [~,new,old] = uniquetol(mesh.X.Values,tol,'ByRows',true,'DataScale', 1) ;
    % Take the barycenter for duplicates groups
        meanMat = sparse(old,1:numel(old),1,numel(new),numel(old)) ;
        meanMat = sparse(1:numel(new),1:numel(new),1./sum(meanMat,2))*meanMat ;
    % Return moved node if needed
        if nargout>1 ; nodesMoved = mesh.X.Values(any(mod(meanMat,1),1),:) ; end
    % Change the node list
        mesh.X.Values = meanMat*mesh.X.Values ;
    % New element list
        mesh.Elems.NodeIdx = mesh.Elems.dataAtIndices(old) ;
    end

end


%% SUBDIVISION, SMOOTHING, ELEMENT ORDER
methods % Functions defined in external files
    % Laplacian Smoothing
    mesh = LaplacianSmoothing(this,lmbda,iterations,tol)
    % Mesh Subdivision
    mesh = CatmullClark(this,iter)
    % Change element types
    mesh = setElementTypes(this,types)
    % Change element order (works for Lagrange elements only)
    mesh = setElementOrder(this,order)
end


%% VIZUALIZATION
methods
    function h = plot(this,varargin)
    % Simple plot of the mesh
        h = pkg.mesh.MeshPlot(varargin{:}) ;
        h.Parent = gca ; 
        h.Mesh = this ;
    end
end




end

