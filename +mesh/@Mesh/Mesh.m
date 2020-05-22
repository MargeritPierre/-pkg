classdef Mesh < matlab.mixin.Copyable
% The general mesh class


%% CONSTRUCTOR / DESTRUCTOR / COPY / SAVE / LOAD
methods
    function this = Mesh(varargin)
    % Class Constructor
        if mod(nargin,2) ; error('wrong number of arguments') ; end
    % Initialize the node position
        this.X = pkg.mesh.MeshFunction('Mesh',this) ;
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
            % Count nodes in each element
                nNodesInElems = this.Elems.nNodes ;
            % Mesh Dimensions
                if isPlanar(this) ; geo = {'tri' 'quad'} ;
                else ; geo = {'tet' 'hex'} ; end
            % Create the list of possible elements
                elemTypes = pkg.mesh.elements.AbstractElement.empty ;
                for gg = 1:numel(geo)
                    for order = 1:3
                        elmt =  pkg.mesh.elements.LagrangeElement(geo{gg},order) ;
                        if ismember(elmt.nNodes,nNodesInElems) ; elemTypes(end+1) = elmt ; end
                    end
                end
            % Assign
                this.Elems.Types = elemTypes ;
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
                    elems = pkg.mesh.elements.ElementTable('Indices',elems,'Types',this.Elems.Types) ;
                else
                    elems = pkg.mesh.elements.ElementTable('NodeIdx',elems,'Types',this.Elems.Types) ;
                end
        end
    % Clean the element list
        elems = clean(elems) ;
    % Set mesh features
        this.Faces = elems.getTableOfUnique('Faces') ;
        this.Edges = elems.getTableOfUnique('Edges') ;
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
end


%% MESH FACES & EDGES
properties (SetAccess = private)
    % Table of faces 
    Faces pkg.mesh.elements.ElementTable
    % Table of edges
    Edges pkg.mesh.elements.ElementTable
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

    function outFace = outerFaces(this)
    % Return an array of logical, true if the face is one an the outer surface
    % The face is on an outer surface if it is linked to only one element
        outFace = sum(this.elem2face,2)==1 ;
    end

    function outNod = outerNodes(this)
    % Return an array of logical, true if the node is on an the outer surface
    % The node is on an outer surface if it is linked to an outer face
        outFace = outerFaces(this) ;
        outNod = logical(this.face2node*outFace) ;
    end

    function bndEdg = boundaryEdges(this)
    % Return an array of logical, true if the edge is on the boundary
    % The edge is on the boundary if it is linked to only one element
        bndEdg = sum(this.elem2edge,2)==1 ;
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
%                         crv{end+1} = nn(1,:) ;
%                         edgBndCrv(bndEdg(1)) = numel(crv) ; 
%                         nn(1,:) = [] ;
%                         bndEdg(1) = [] ;
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

    function [normals,frames] = faceNormals(this)
    % Return an array containing the face normals
    % NORMALS size: [nFaces 3]
    % Also return the face frames if needed
    % FRAMES size: [nFaces 3 3] ([iface ivect icoord])
        % First check if the mesh is planar
            [ispl,normal,~,frame] = isPlanar(this) ;
            if ispl 
                normals = repmat(normal(:)',[this.nFaces 1]) ;
                frames = repmat(reshape(frame,[1 3 3]),[this.nFaces 1 1]) ;
                return
            end
        % Otherwise fit a plane on each face
            normals = zeros(this.nFaces,3) ;
            frames = zeros(this.nFaces,3,3) ;
            C = this.circumCenters(this.Faces) ;
            for ff = 1:this.nFaces
                ii = this.Faces(ff) ;
                ii = ii(~isnan(ii)) ;
                P = this.Nodes(ii,:) - C(ff) ;
                % Compute the coordinates covariance
                    [V,s] = eig(P'*P,'vector') ;
                % V is the frame basis, s the variances
                    [~,ind] = sort(s,'asccend') ;
                    V = V(:,ind) ;
                % Return the objects
                    frames(ff,:,:) = V' ;
                    normals(ff,:) = frames(ff,end,:) ;
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
            isP1Quad = mesh.Elems.nNodes == 4 ;
            mesh.Elems.NodeIdx(isP1Quad,:) = mesh.Elems.NodeIdx(isP1Quad,[1 2 4 3]) ;
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
                notConverged = NaN ;
                while ~isempty(notConverged)
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

    function mesh = replicate(this,N)
    % Reproduce the mesh N times
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        if N<2 ; return ; end
        mesh.X = repmat(this.X.Values,[N 1]) ;
        elems = repmat(this.Elems,[N 1])  ; 
        elems.NodeIdx = elems.NodeIdx + kron((0:N-1)',ones(this.nElems,1)) ;
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
                F = [F O ; permute(O,[2 1 3]) O(1:dC,1:dC,:)] ;
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
            F = [cos(THETA) sin(THETA) ; -sin(THETA) cos(THETA)] ;
            if mesh.nCoord>2 ; F(end+1,end+1,:) = 1 ; end
        % Rotate
            mesh.applyTransform(F) ;
        % Re-move the mesh
            mesh.move(CEN) ;
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
            remove = full(sparse(remove(:),ones(numel(remove),1),true,mesh.Elems,1)) ; 
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
        mesh.Elems.NodeIdx = mesh.Elems.dataAtIndices(old);
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
    function h = plot(this)
    % Simple plot of the mesh
        h = pkg.mesh.MeshPlot() ;
        h.Parent = gca ; 
        h.Mesh = this ;
    end
end




end

