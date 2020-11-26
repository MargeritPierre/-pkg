classdef Mesh < matlab.mixin.Copyable ...
                    & pkg.geometry.NodeObject ...
                    & pkg.geometry.mesh.elements.TableObject
% The general mesh class


%% CONSTRUCTOR / DESTRUCTOR / COPY / SAVE / LOAD
methods
    function this = Mesh(varargin)
    % Class Constructor
    % Process inputs
        if nargin==1 % mesh from single input (see below)
            this = this.meshFromOneInput(varargin{:}) ; 
            return ; 
        end
        if mod(nargin,2) ; error('wrong number of arguments') ; end
    % Input arguments
        PropNames = varargin(1:2:end-1) ;
        PropValues = varargin(2:2:end) ;
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
                    nCoord = size(x,2) ;
                    nodeIdx = [1:size(x,1)-1 ; 2:size(x,1)]' ;
                    types = pkg.geometry.mesh.elements.base.Bar ;
                case 3 % surface mesh
                    [nY,nX,nCoord] = size(x) ;
                    p1 = (1:nY-1)'+ nY*(0:nX-2) ; p2 = p1+nY ; p3 = p2+1 ; p4 = p1+1 ;
                    nodeIdx = [p1(:) p2(:) p3(:) p4(:)] ;
                    types = pkg.geometry.mesh.elements.base.Quadrangle ;
                case 4 % volume mesh
                    error('Automatic volume mesh not implemented yet') ;
            end
            this.Nodes = reshape(x,[],nCoord) ;
            idx = padarray(nodeIdx,[0 1],1,'pre') ;
            this.Elems = pkg.geometry.mesh.elements.ElementTable('Types',types,'Indices',idx) ;
        else
    % Mesh from a specific object class
            switch class(input)
                case 'pkg.geometry.mesh.Mesh'
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
    end
end

%% ELEMENT ASSIGNMENT
methods
    function assignElements(this,types)
    % Automatically assign elements types to the mesh
        if this.nElems==0 ; return ; end
        if nargin<2 ; types = 'auto' ; end
    % Process types options
        if ischar(types)
            switch types
                case '1D'
                case '2D'
                case '3D'
                case 'auto' % Preset list of base elements
                    if isPlanar(this) ; types = '2D' ;
                    else ; types = '3D' ;
                    end
                otherwise
                    error('Unsupported option for element types.') ;
            end
        end
    % Apply
        assignElements@pkg.geometry.mesh.elements.TableObject(this,types) ;
    end
end


%% OPERATIONS ON CONNECTIVITIES (node/edge/face splitting, merging, etc)
methods
    function mesh = splitNodes(this,nod)
    % Split the mesh at given nodes
    % Create a new node for each attached element
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        if islogical(nod) ; nod = find(nod) ; end
        nod = unique(nod(:)) ;
        for nn = nod(:)'
            ind = ismember(mesh.Elems.NodeIdx,nn) ;
            newNodes = [nn mesh.nNodes + (1:sum(ind(:))-1)] ;
            mesh.Elems.NodeIdx(ind(:)) = newNodes(:) ;
            mesh.Nodes(newNodes,:) = repmat(mesh.Nodes(nn,:),[numel(newNodes) 1]) ;
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
        x = this.Nodes ;
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
        if nargout==0 ; this.Nodes = x ; this.Elems = elems ; 
        else ; mesh = pkg.geometry.mesh.Mesh('Nodes',x,'Elems',elems) ; 
        end
    end
end


%% GEOMETRY
methods
    function bool = isSurface(this)
    % Is this mesh representing a surface ?
        bool = all([this.Elems.Types.nDims]==2) ;
    end
    
    function bool = isVolume(this)
    % Is this mesh representing a volume ?
        bool = all([this.Elems.Types.nDims]==3) ;
    end
    
    function C = centroid(this,table)
    % Return the centroids associated to the given list
    % By default, list is this.Elems
        if nargin<2 ; table = this.Elems ; end
        if ischar(table) ; table = this.(table) ; end
        C = table.meanDataAtIndices(this.Nodes) ;
        C = reshape(C,table.nElems,this.nCoord) ;
    end
    
    function [Ne,ed,Nn,no] = boundaryNormals(this)
    % Return the normals at the boundaries 
        [ed,~,no] = boundaryFeatures(this) ;
        Ne = this.getNormals(this.Edges.subpart(ed)) ;
        if nargout>=3 % boundary node normals: mean of attached boundary edges
            e2n = this.edge2node ; % connectivity
            e2n = e2n(no,ed) ; % only between boundary features
            Nn = (e2n*Ne)./sum(e2n,2) ; % mean over attached edges
            Nn = Nn./sqrt(sum(Nn.^2,2)) ;
        end
    end
    
    function normals = nodeNormals(this)
    % (Mandatory for pkg.geometry.NodeObject)
    % Return the 3D normal to the mesh evaluated on nodes
    % normals [nNodes nCoord]
    % From element normals to node normals
        normals = this.elem2node('mean')*this.getNormals(this.Elems) ;
    % Normalize
        normals = normals./sqrt(sum(normals.^2,2)) ;
    end
    
    function normals = getNormals(this,features,E)
    % Return the local normals of features evaluated at local coordinates E
    % features: pkg.geometry.mesh.elements.ElementTable [nFeat nMaxNodesByFeat]
    % E: [nFeat nMaxDims]
    % normals: [nFeat 3]
        if nargin<2 ; features = this.Edges ; end
        if ischar(features) ; features = this.(features) ; end
        if nargin<3 ; E = 'centroid' ; end
        frames = permute(this.getFrames(features,E),[3 2 1]) ;
    % Retrieve the normal
        normals = NaN(features.nElems,3) ;
        for tt = 1:numel(features.Types)
            elmtIdx = features.TypeIdx==tt ;
            elmtType = features.Types(tt) ;
            nDims = elmtType.nDims ;
            if nDims==1
                normals(elmtIdx,:) = -frames(elmtIdx,:,2) ;
            else
                normals(elmtIdx,:) = frames(elmtIdx,:,3) ;
            end
        end
    end
    
    function frames = getFrames(this,features,E)
    % Return the local frames of features evaluated at local coordinates E
    % features: pkg.geometry.mesh.elements.ElementTable [nFeat nMaxNodesByFeat]
    % E: [nFeat nMaxDims]
    % frames: [3 3 nFeat] (ivect icoord ifeat)
        if nargin<2 ; features = this.Elems ; end
        if ischar(features) ; features = this.(features) ; end
        if nargin<3 ; E = 'centroid' ; 
        elseif isnumeric(E) && size(E,1)~=features.nElems ; error('Wrong input format') ; end
        frames = NaN(3,3,features.nElems) ;
    % Node coordinates for each element
        Xf = features.dataAtIndices(this.Nodes) ;
        Xf = padarray(Xf(:,:,1:min(end,3)),[0 0 3-size(Xf,3)],0,'post') ; % force 3D coordinates
    % Dispatch to element types
        for tt = 1:numel(features.Types)
            elmtType = features.Types(tt) ; 
            elmtIdx = features.TypeIdx==tt ;
            nE = sum(elmtIdx) ;
        % Local coordinates
            if strcmp(E,'centroid') ; ee = repmat(elmtType.centroid,[nE 1]) ;
            else ; ee = E(elmtIdx,:) ; end
        % Jacobian & gradient
            dN_de = elmtType.evalJacobianAt(ee) ; % [nE nNodes nDims]
            dX_de = sum(dN_de.*permute(Xf(elmtIdx,1:elmtType.nNodes,:),[1 2 4 3]),2) ; % [nE 1 nDims nCoord]
            dX_de = permute(dX_de,[4 3 1 2]) ; % [nCoord nDims nE]
        % Local frame
            frames(1:elmtType.nDims,:,elmtIdx) = permute(dX_de,[2 1 3]) ; 
        end
    % Fill with other vectors if needed
        % 1D elements
            is1D = sum(~isnan(frames(:,1,:)),1)==1 ;
            isV1z = all(abs(frames(1,:,is1D))==[0 0 1],2) ; % is the tangent along Z ?
            frames(3,:,is1D) = [0 0 1].*~isV1z + [0 1 0].*isV1z ; % third vector = Z or Y by default
            frames(2,:,is1D) = pkg.math.vectprod(frames(3,:,is1D),frames(1,:,is1D),2) ; % second vector is the vectorial product
        % 2D elements
            is2D = sum(~isnan(frames(:,1,:)),1)==2 ;
            frames(3,:,is2D) = pkg.math.vectprod(frames(1,:,is2D),frames(2,:,is2D),2) ; % third vector is the vectorial product
    % Normalize
        frames = frames./sqrt(sum(frames.^2,2)) ;
    end

    function mesh = sortElems(this,dir)
    % Sort the element nodes of a planar mesh in clockwise 
    % or counter-clockwise(==trigo) direction
    % WARNING: only works with 1st order TRIS or QUADS !!!
        if ~this.isPlanar ; error('Cannot sort nodes of a non-planar mesh') ; end
        if nargin<2 ; dir = 'trigo' ; end
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        % Barycentered coordinates
            x = mesh.Elems.dataAtIndices(mesh.Nodes) ;
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
    
    function sz = elemSize(this,features)
    % Compute the typical size of each element
    % 1D: length , 2D: area , 3D: volume
    % <TODO> Generalize this function
        if nargin<2 ; features = this.Elems ; end
        [E,ie] = features.getListOf('GaussIntegrationPoints') ;
        J = this.detJacobian(E,ie,features) ;
        [W,~] = features.getListOf('GaussIntegrationWeights') ;
        sz = abs(J).*W ;
        sz = accumarray(ie(:),sz(:),[features.nElems 1]) ;
    end
    
end


%% POINT LOCALIZATION
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
        if nargin<4 ; X = this.Nodes ; end
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
    % Points in element bounding box ?
        [pp,ee] = inElmtBBox(this,P,features,X,tol) ;
        if bboxOnly ; M = logical(sparse(pp,ee,1,nPts,features.nElems)) ; return ; end
    % Localize each point in each feature
        [uE,~,dist] = localize(this,P,features,false,X,tol,pp,ee) ;
        inside = ~any(isnan(uE),2) & dist<=tol  ;
        pp = pp(inside) ; ee = ee(inside) ;
    % Build the matrix
        M = logical(sparse(pp,ee,1,nPts,features.nElems)) ;
    end
    
    function [pp,ee] = inElmtBBox(this,P,features,X,tol,ip,ie)
    % Return indices pp and ee corresponding to P(pp) inside elmt(ee)
    % P = [nPts nCoord]
    % X = [nNodes nCoord]
    % By default, each point is tested against each element of the mesh
    % Optionally, (ip ie) [nTest 1] can be provided to test specific point with
    % specific elements (see below)
        if nargin<3 ; features = this.Elems ; end
        if nargin<4 ; X = this.Nodes ; end
        if nargin<5 ; tol = this.defaultTolerance(X) ; end
        nPts = size(P,1) ;
    % Element node coordinates
        Xe = features.dataAtIndices(X) ; % [nElems nMaxNodesByElem nCoord]
        if nargin>5 % test specific points vs specific elements
        % Put the coordinates to the third dimension
            P = reshape(P,nPts,1,[]) ; % [nPts 1 nCoord]
            if nargin>5 ; P = P(ip,:,:) ; end % [nTest 1 nCoord]
        % Compute bounding boxes
            Xe = Xe(ie,:,:) ; % [nTest nMaxNodesByElem nCoord] 
            elmtBBox = [min(Xe,[],2) max(Xe,[],2)] ; % [nTest 2 nCoord]
        % Test if point inside
            in = all(P>=elmtBBox(:,1,:)-tol & P<=elmtBBox(:,2,:)+tol,3) ; % [nTest 1]
            ii = find(in) ;
            pp = ip(ii) ; ee = ie(ii) ;
        else % test all points vs all elements...
            elmtBBox = [min(Xe,[],2) max(Xe,[],2)] ; % [nElems 2 nCoord]
            [pp,ee] = pkg.data.inDomain(P,elmtBBox) ; 
        end
    end
    
    function [E,ie,dist] = localize(this,P,features,extrap,X,tol,ip,ie)
    % Localize points P on the mesh. 
    % Returns the local coordinates corresponding to the feature's closest point
    % Find E = argmin(norm(features.evalAt(E)*X-P))
    % P = [nPts nCoord]
    % features = ElementTable: mesh Elems (default), Faces or Edges
    % extrap = return point localization even outside the features (else E(outside,:) = NaN)
    % extrap is by default false as (ip,ie) are not provided
    % X = [nNodes nCoord] (mesh.X by default)
    % tol: convergence tolerance
    % ip,ie = [nLocal 1] optionnal: restricts the test between P(ip) and features(ie)
    % E = [nPts nElems nMaxElemDim] or [nLocal nMaxElemDim]
        if nargin<3 ; features = this.Elems ; end
        if nargin<4 ; extrap = false ; end
        if nargin<5 ; X = this.Nodes ; end
        if nargin<6 ; tol = this.defaultTolerance(X) ; end
        if nargin==7 ; error('Point AND element indices couples (ip,ie) must be provided') ; end
        itMax = 10 ;
        debug = false ;
        if debug ; pl0 = plot(P(:,1),P(:,2),'or') ; pl = plot(P(:,1)*NaN,P(:,2)*NaN,'.b') ; end
    % Localization mode
        if nargin>7 % Do what is asked, all inputs are available
            mode = 'custom' ;
        elseif extrap && nargin<7 % Localize each point in each mesh feature
            mode = 'closestFeatures' ;
        else % Localize each point in the closest mesh features (bbox only)
            mode = 'inBBox' ;
        end
    % If extrapolation is allowed, localize each point in each mesh feature
        nPts = size(P,1) ; nElems = features.nElems ;
        switch mode
            case 'custom' % Do what is asked, all inputs are available
            case 'allInAll' % Localize each point in each feature
                ip = repmat(1:nPts,[1 nElems]) ; 
                ie = repelem(1:nElems,nPts*ones(1,nElems)) ; 
            case 'closestFeatures' % Localize each point only in the closest features
                nn = this.closestNode(P,X) ;
                f2n = features.sparse ;
                [ip,ie] = find(f2n(nn,:)) ;
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
        dist = NaN(nLocal,1) ;
    % For each element type
        for typeIdx = 1:features.nTypes
            elmtType = features.Types(typeIdx) ; % The type of feature to deal with
            % Find the corresponding element indices
                elmtIdx = find(features.TypeIdx==typeIdx) ; % The indices in the feature list
                [isType,ii] = ismember(ie,elmtIdx) ; % Keep only elements of interest
                isType = find(isType) ;
                elmtIdx = elmtIdx(ii(isType)) ; % The element indices we will have to deal with
                nE = length(elmtIdx) ;
            % Retrieve the element node coordinates
                xe = Xe(elmtIdx,1:elmtType.nNodes,:) ; % [nE nElmtNodes nCoord]
                xe = permute(xe,[3 1 2]) ; % [nCoord nE nElmtNodes]
            % Intialize the algorithm
                e0 = mean(elmtType.NodeLocalCoordinates,1) ;
                if norm(reshape(elmtType.evalHessianAt(e0),[],1))<tol ; methodOrder = 1 ;
                else ; methodOrder = 2 ; end
                e = repmat(e0,[nE 1]) ; % [nE nElmtNodes]
            % Loop until norm(x-P)<tol
                res = NaN(this.nCoord,nE) ;
                notConverged = 1:nE ; it = 0 ;
                while ~isempty(notConverged) && it<itMax
                    it = it + 1 ;
                    e_nc = e(notConverged,:) ; % [nE_nc nElmtDims]
                    xe_nc = xe(:,notConverged,:) ;
                % Shape function evalutation (+derivatives)
                    N = elmtType.evalAt(e_nc) ; % [nE_nc nElmtNodes]
                    dN_de = elmtType.evalJacobianAt(e_nc) ; % [nE_nc nElmtNodes nElmtDims]
                    if methodOrder>1 ; d2N_de2 = elmtType.evalHessianAt(e_nc) ; end % [nE_nc nElmtNodes nElmtDims nElmtDims]
                % Function guess evaluation (+derivatives)
                    x = sum(permute(N,[3 1 2]).*xe_nc,3) ; % [nCoord nE_nc]
                    if debug ; pl.XData = x(1,:) ; pl.YData = x(2,:) ; drawnow ; end
                    dx_de = sum(permute(dN_de,[4 1 2 3]).*xe_nc,3) ; % [nCoord nE_nc 1 nElmtDims]
                    if methodOrder>1 ; d2x_de2 = sum(permute(d2N_de2,[5 1 2 3 4]).*xe_nc,3) ; end % [nCoord nE_nc 1 nElmtDims nElmtDims]
                    dx_de = permute(dx_de,[1 4 2 3]) ; % [nCoord nElmtDims nE_nc]
                    if methodOrder>1 ; d2x_de2 = permute(d2x_de2,[1 4 5 2 3]) ; end % [nCoord nElmtDims nElmtDims nE_nc]
                % Residue check
                    r = x-P(:,isType(notConverged)) ; % [nCoord nE_nc]
                    res(:,notConverged) = r ; % [nCoord nE]
                    stillNotConverged = sum(r.^2,1)>tol^2 ;
                    r = r(:,stillNotConverged) ;
                    notConverged = notConverged(stillNotConverged) ;
                    if isempty(notConverged) ; break ; end
                % Update points where it needs to be
                    J = dx_de(:,:,stillNotConverged) ; % [nCoord nElmtDims nE_nc]
                    if methodOrder==1
                        de = pkg.math.mldivide(J,r) ;
                    else
                        J2 = pkg.math.mtimes(permute(J,[2 1 3]),J) ; % [nElmtDims nElmtDims nE_nc]
                        Jr = pkg.math.mtimes(permute(J,[2 1 3]),permute(r,[1 3 2])) ; % [nElmtDims 1 nE_nc]
                        H = sum(permute(r,[1 3 4 2]).*d2x_de2(:,:,:,stillNotConverged),1) ; % [1 nElmtDims nElmtDims nE_nc] ;
                        H = J2 + permute(H,[2 3 4 1]) ; % [nElmtDims nElmtDims nE_nc] ;
                        de = pkg.math.mldivide(H,Jr) ;
                    end
                    e(notConverged,:) = e(notConverged,:) - de(:,:).' ;
                end
            % Assign
                E(isType,:) = e ;
                dist(isType) = sqrt(sum(res.^2,1)) ;
            % If not extrapolating, check that E correspond to local
            % coordinate INSIDE the element
                if ~extrap
                    outside = ~elmtType.isInside(E(isType,:)) ;
                    E(isType(outside),:) = NaN ;
                    dist(isType(outside)) = NaN ;
                end
        end
    % Keep only one localization point
        switch mode
            case 'allInAll' % each point is localized in each element
                E = reshape(E,nPts,nElems,[]) ; 
                ie = reshape(ie,nPts,nElems,[]) ; 
                dist = reshape(dist,nPts,nElems,[]) ; 
            case 'custom' % let it as it has been asked, do nothing
            otherwise % keep only the closest element
                [~,ind] = sort(dist,'ascend') ;
                [ip,ia] = unique(ip(ind),'first') ;
                ind = ind(ia) ;
                uE = NaN(nPts,size(E,2)) ; uE(ip,:) = E(ind,:) ; E = uE ;
                uie = ones(nPts,1) ; uie(ip) = ie(ind) ; ie = uie ;
                udist = NaN(nPts,1) ; udist(ip) = dist(ind) ; dist = udist ;
        end
        if strcmp(mode,'allInAll') ; E = reshape(E,nPts,nElems,[]) ; end
        if debug ; delete(pl0) ; delete(pl) ; end
    end
end

%% MESH INTERPOLATION & DIFFERENCIATION
methods
% FOR THE FOLLOWING FUNCTIONS:
% f_n are the value of any function at nodes
% if ie is NOT given (or empty), P_or_E is global and is localized in the mesh first
% if ie IS given, P_or_E is condidered to be local coordinates
    
    function M = interpMat(this,P_or_E,ie,features,extrap,X,tol)
    % Return a sparse interpolation matrix so that f(P) = M*f_n
        if nargin<2 ; P_or_E = this.centroid ; end
        if nargin<3 ; ie = [] ; end
        if nargin<4 ; features = this.Elems ; end
        if nargin<5 ; extrap = false ; end
        if nargin<6 ; X = this.Nodes ; end
        if nargin<7 ; tol = this.defaultTolerance(X) ; end
        M = sparse(size(P_or_E,1),this.nNodes) ;
    % Localize the given points in the mesh if needed
        if isempty(ie)  ; [E,ie] = this.localize(P_or_E,features,extrap,X,tol) ;
        else ; E = P_or_E ; end
        features = features.subpart(ie) ; 
    % Evaluate the element shape functions at the given coordinates
        for typeIdx = 1:features.nTypes
            elmtType = features.Types(typeIdx) ;
            elmtIdx = find(features.TypeIdx==typeIdx) ;
            N = elmtType.evalAt(E(elmtIdx,1:elmtType.nDims)) ;
            M = M + sparse(repmat(elmtIdx(:),[1 elmtType.nNodes]),double(features.NodeIdx(elmtIdx,1:elmtType.nNodes)),N,size(P_or_E,1),this.nNodes) ;
        end
    end
    
    function M = gradMat(this,P_or_E,ie,features,extrap,X,tol)
    % Return sparse matrices M so that df(P)/dx_i = M{i}(P)*f_n
        if nargin<3 ; ie = [] ; end
        if nargin<4 ; features = this.Elems ; end
        if nargin<5 ; extrap = false ; end
        if nargin<6 ; X = this.Nodes ; end
        if nargin<7 ; tol = this.defaultTolerance(X) ; end
    % Localize the given points in the mesh
        if nargin<2 % element centroids
            E = features.getListOf('centroid') ;
        else % given points
            if isempty(ie)  ; [E,ie] = this.localize(P_or_E,features,extrap,X,tol) ;
            else ; E = P_or_E ; end
            features = features.subpart(ie) ; 
        end
    % Evaluate the element shape function gradient at the given coordinates
        M = repmat({sparse(size(E,1),this.nNodes)},[this.nCoord 1]) ;
        xe = features.dataAtIndices(X) ; % [nE nMaxNodesByElmt nCoord] 
        for typeIdx = 1:features.nTypes
            elmtType = features.Types(typeIdx) ;
            elmtIdx = find(features.TypeIdx==typeIdx) ;
            dN_de = elmtType.evalJacobianAt(E(elmtIdx,1:elmtType.nDims)) ; % [nE nElmtNodes nElmtDims]
            dx_de = permute(sum(xe(elmtIdx,1:elmtType.nNodes,:).*permute(dN_de,[1 2 4 3]),2),[1 3 4 2]) ; % [nE nCoord nElmtDims] 
            de_dx = pkg.math.pinv(permute(dx_de,[2 3 1])) ; % [nElmtDims nCoord nE]
            dN_dx = permute(pkg.math.mtimes(permute(dN_de,[2 3 1]),de_dx),[3 1 2]) ; % [nE nElmtNodes nCoord]
            iii = repmat(elmtIdx(:),[1 elmtType.nNodes]) ;
            jjj = double(features.NodeIdx(elmtIdx,1:elmtType.nNodes)) ;
            for cc = 1:this.nCoord
                M{cc} = M{cc} + sparse(iii,jjj,dN_dx(:,:,cc),size(E,1),this.nNodes) ;
            end
        end
    end
    
    function M = grad2Mat(this,P_or_E,ie,features,extrap,X,tol)
    % Return sparse matrices M so that d2f(P)/(dx_i.dx_j) = M{i,j}(P)*f_n
    % d2f_dx2 = (d_dx)(df_de.de_dx) = d2f_de2.de_dx.de_dx + df_de.d2e_dx2
    % de_dx is the inverse Jacobian
    % d2e_dx2 is the inverse Hessian
        if nargin<3 ; ie = [] ; end
        if nargin<4 ; features = this.Elems ; end
        if nargin<5 ; extrap = false ; end
        if nargin<6 ; X = this.Nodes ; end
        if nargin<7 ; tol = this.defaultTolerance(X) ; end
    % Localize the given points in the mesh
        if nargin<2 % element centroids
            E = features.getListOf('centroid') ;
        else % given points
            if isempty(ie)  ; [E,ie] = this.localize(P_or_E,features,extrap,X,tol) ;
            else ; E = P_or_E ; end
            features = features.subpart(ie) ; 
        end
    % Evaluate the element shape function gradient at the given coordinates
        M = repmat({sparse(size(E,1),this.nNodes)},[this.nCoord this.nCoord]) ;
        xe = features.dataAtIndices(X) ; % [nE nMaxNodesByElmt nCoord] 
        for typeIdx = 1:features.nTypes
            elmtType = features.Types(typeIdx) ;
            elmtIdx = find(features.TypeIdx==typeIdx) ;
            % Iverse Jacobian
                dN_de = elmtType.evalJacobianAt(E(elmtIdx,1:elmtType.nDims)) ; % [nE nElmtNodes nElmtDims]
                dx_de = permute(sum(xe(elmtIdx,1:elmtType.nNodes,:).*permute(dN_de,[1 2 4 3]),2),[1 3 4 2]) ; % [nE nCoord nElmtDims] 
                de_dx = pkg.math.pinv(permute(dx_de,[2 3 1])) ; % [nElmtDims nCoord nE]
            % Hessian
                d2N_de2 = elmtType.evalHessianAt(E(elmtIdx,1:elmtType.nDims)) ; % [nE nElmtNodes nElmtDims nElmtDims]
                %d2x_de2 = permute(sum(xe(elmtIdx,1:elmtType.nNodes,:).*permute(d2N_de2,[1 2 5 3 4]),2),[1 3 4 5 2]) ; % [nE nCoord nElmtDims nElmtDims]
            % Second-order derivatives
                d2N_de2 = permute(d2N_de2,[3 4 1 2]) ; % [nElmtDims nElmtDims nE nElmtNodes]
                d2N_dxde = pkg.math.mtimes(d2N_de2,de_dx) ; % [nCoord nCoord nE nElmtNodes]
                dN2_dx2 = pkg.math.mtimes(permute(de_dx,[2 1 3]),d2N_dxde) ; % [nCoord nCoord nE nElmtNodes]
                dN2_dx2 = permute(dN2_dx2,[3 4 1 2]) ; % [nE nElmtNodes nCoord nCoord]
            % Sparse matrices
                iii = repmat(elmtIdx(:),[1 elmtType.nNodes]) ;
                jjj = double(features.NodeIdx(elmtIdx,1:elmtType.nNodes)) ;
                for c1 = 1:this.nCoord
                    for c2 = 1:this.nCoord
                        M{c1,c2} = M{c1,c2} + sparse(iii,jjj,dN2_dx2(:,:,c1,c2),size(E,1),this.nNodes) ;
                    end
                end
        end
    end
    
    function J = detJacobian(this,P,ie,features,extrap,X,tol)
    % Return a sparse matrix M so that int(f) = M*f_n
    % with f_n the value of any function at nodes
    % if ie is NOT given (or empty), P is localized in the mesh first
    % if ie IS given, P is condidered to be local coordinates (E)
        if nargin<3 ; ie = [] ; end
        if nargin<4 ; features = this.Elems ; end
        if nargin<5 ; extrap = false ; end
        if nargin<6 ; X = this.Nodes ; end
        if nargin<7 ; tol = this.defaultTolerance(X) ; end
    % Localize the given points in the mesh
        if nargin<2 % element centroids
            [E,ie] = features.getListOf('GaussIntegrationPoints') ;
        else % given points
            if isempty(ie)  ; [E,ie] = this.localize(P,features,extrap,X,tol) ;
            else ; E = P ; end
        end
        features = features.subpart(ie) ; 
    % Compute the jacobian and its determinant
        Xe = features.dataAtIndices(X) ;
        J = zeros(features.nElems,1) ;
        for tt = 1:features.nTypes
            elmtIdx = features.TypeIdx==tt ;
            elmtType = features.Types(tt) ;
            dN_de = elmtType.evalJacobianAt(E(elmtIdx,:)) ; % [nE nNodesInElmt nElmtDims]
            dX_de = sum(Xe(elmtIdx,1:elmtType.nNodes,:).*permute(dN_de,[1 2 4 3]),2) ; % [nE 1 nCoord nElmtDims]
            if elmtType.nDims==this.nCoord
                J(elmtIdx) = pkg.math.det(permute(dX_de,[3 4 1 2])) ; % [1 1 nE 1] ;
            else % (for nCoord~=nElmtDims) use the squared jacobian 
                dX_de2 = sum(dX_de.*permute(dX_de,[1 4 3 2]),3) ; % [nE nElmtDims 1 nElmtDims]
                J(elmtIdx) = sqrt(pkg.math.det(permute(dX_de2,[2 4 1 3]))) ; % [1 1 nE 1] ;
            end
        end
    end
    
    function [P,ie,E] = gaussIntegrationPoints(this,features,X)
    % Return all gauss points of the features
    % P: [nF nCoord], ie: [nF 1], E: [nF maxFeatureDim]
        if nargin<2 ; features = this.Elems ; end
        if nargin<3 ; X = this.Nodes ; end
        [E,ie] = features.getListOf('GaussIntegrationPoints') ;
        M = interpMat(this,E,ie,features) ;
        P = M*X ;
    end
    
end

%% MESH OPERATIONS (MOTION, BOOLEAN)
methods 

    function mesh = plus(A,B)
    % Addition of two meshes (NEED CLEANUP TO WELD THE MESHES!)
        mesh = copy(A) ; % Always take a copy
        % New mesh coordinates values
            nNodes = A.nNodes + B.nNodes ;
            x = NaN(nNodes,max(A.nCoord,B.nCoord)) ;
            x(1:A.nNodes,1:A.nCoord) = A.Nodes ;
            x(A.nNodes+1:end,1:B.nCoord) = B.Nodes ;
        % New element table
            elems = [A.Elems B.Elems] ;
            elems.NodeIdx(A.nElems+1:end,:) = elems.NodeIdx(A.nElems+1:end,:) ...
                                            + uint32(A.nNodes.*logical(elems.NodeIdx(A.nElems+1:end,:))) ;
        % Return/Modify the mesh
            mesh.Nodes = x ; 
            mesh.Elems = elems ;
    end
    
    function meshes = splitParts(this)
    % Split the independent parts of a mesh into several-meshes
        e2n = this.elem2node ;
        e2e = e2n'*e2n ;
        meshes = pkg.geometry.mesh.Mesh.empty ;
        elems = this.Elems ;
        while ~isempty(e2e)
            belong = e2e\sparse(1,1,1,size(e2e,1),1)~=0 ;
            belongElems = pkg.geometry.mesh.elements.ElementTable('Types',this.Elems.Types,'Indices',this.Elems.Indices(belong,:)) ;
            meshes(end+1) = pkg.geometry.mesh.Mesh('Nodes',this.Nodes,'Elems',belongElems) ;
            elems.Indices = elems.Indices(~belong,:) ;
            e2e = e2e(~belong,~belong) ;
        end
    end
    
    function mesh = merge(A,B,tol)
    % Merge two meshes (Mandatory for pkg.geometry.NodeObject)
        if nargin<3 ; tol = A.defaultTolerance ; end
    % Merge nodes
        [mesh,nodA,nodB] = merge@pkg.geometry.NodeObject(A,B,tol) ;
    % Merge elements
        elems = [A.Elems B.Elems] ;
        nodA = [0 ; nodA(:)] ; nodB = [0 ; nodB(:)] ;
        elems.NodeIdx(1:A.nElems,:) = nodA(elems.NodeIdx(1:A.nElems,:)+1) ;
        elems.NodeIdx(A.nElems+(1:B.nElems),:) = nodB(elems.NodeIdx(A.nElems+(1:B.nElems),:)+1) ;
        mesh.Elems = elems ;
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
        if numel(x)==1 ; x = repmat(this.Nodes,[1 1 x]) ; end
        if size(x,1)~=size(this.Nodes,1) ; error('Wrong format for input X.') ; end
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
        indices = permute(indices,[1 3 2]) ; % [this.nElems nRep nMaxNodes+1]
    % Fill ?
        types = this.Elems.Types ;
        if fill % need to return extruded elements
            if loop ; indRep = [1:nRep ; [2:nRep,1]] ;
            else ; indRep = [1:nRep-1 ; 2:nRep] ;
            end
            nElems = this.nElems*size(indRep,2) ; % final number of elements
            nMaxNodes = this.nMaxNodesByElem ; % current maximum number of nodes
            indices = [ reshape(indices(:,1:end-1,1),nElems,1) ...
                        reshape(indices(:,indRep(1,:),2:end),nElems,nMaxNodes) ...
                        reshape(indices(:,indRep(2,:),2:end),nElems,nMaxNodes) ...
                        ] ;
            for tt = 1:numel(types)
            % Valid nodes (for heterogeeous element types...
                validNodes = (1:types(tt).nNodes)'+[0 1]*nMaxNodes ;
            % Depending on the initial element type...
                switch class(types(tt))
                    case 'pkg.geometry.mesh.elements.base.Bar'
                        types(tt) = pkg.geometry.mesh.elements.base.Quadrangle ;
                        validNodes(:,2) = flip(validNodes(:,2)) ; % quad numerotation...
                    case 'pkg.geometry.mesh.elements.base.Triangle'
                        types(tt) = pkg.geometry.mesh.elements.base.Prism ;
                    case 'pkg.geometry.mesh.elements.base.Quadrangle'
                        types(tt) = pkg.geometry.mesh.elements.base.Hexahedron ;
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
            indices = reshape(indices,this.nElems*nRep,[]) ; % to 2d array of indices
        end
        elems = pkg.geometry.mesh.elements.ElementTable('Types',types,'Indices',indices) ;
    % Assign
        if nargout==0 ; this.Nodes = x ; this.Elems = elems ; 
        else ; mesh = pkg.geometry.mesh.Mesh('Nodes',x,'Elems',elems) ; 
        end
    end  
    
end


%% MESH SLICING, CUTTING, SPLITTING,... BY A LEVEL SET
methods
    
    function sBool = lvlSetSign(this,fcnX,tol,X)
    % Evaluate a levelset f = fcn(X) on mesh nodes X and return a signed
    % boolean: -1 (inside), 0 (on), 1 (outside)
        if nargin<3 ; tol = this.defaultTolerance ; end
        if nargin<4 ; X = this.Nodes ; end
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
        if nargin<5 ; X = this.Nodes ; end
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
        if nargin<5 ; X = this.Nodes ; end
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
    
    function mesh = cut(this,fcn,tol,X)
    % Return a cut of the mesh so that for any X in mesh, fcn(X)<=0
        if nargin<3 ; tol = this.defaultTolerance ; end
        if nargin<4 ; X = this.Nodes ; end
    % Evaluate the levelset on nodes
        fcnX = fcn(X) ;
    % Get which elements are crossing the level set
        [fBool,eBool] = featureLevelSetSign(this,fcnX,this.Elems,tol,X) ;
        indElemsOnLvlSt = find(fBool==0) ;
        eBool = eBool(indElemsOnLvlSt,:) ;
        elemsOnLvlSt = this.Elems.subpart(indElemsOnLvlSt) ;
    % Edges on the level set: we need to retrieve the connectivity...
        [edg,elmt,elmtIdx] = find(this.ElemEdges(:,indElemsOnLvlSt)) ;
        edgesOnLvlSt = full(sparse(elmt,elmtIdx,edg)) ; % [nElemsOnLvlSt maxEdgByElem]
    % Keep useful information
        types = elemsOnLvlSt.Types ;
        typeIdx = elemsOnLvlSt.TypeIdx ;
        nodeIdx = [elemsOnLvlSt.NodeIdx edgesOnLvlSt+this.nNodes] ; % indices of nodes & edges
        nMaxNodesByElem = size(elemsOnLvlSt.NodeIdx,2) ;
    % Slice elements crossing the levelset
        sliceElems = struct('IN',[],'ON',[],'OUT',[]) ;
        for tt = 1:elemsOnLvlSt.nTypes
        % Retrieve concerned elements
            elmtIdx = find(typeIdx==tt) ;
            if isempty(elmtIdx) ; continue ; end
        % Get the slicing cases structure
        % see pkg.geometry.mesh.elements.base.BaseElement.slicecases
            S = types(tt).sliceCases ;
        % Get which element correspond to which slicing configuration
            tests = cat(1,S.Test) ;
            [~,sliceConfig] = ismember(eBool(elmtIdx,1:types(tt).nNodes),tests(:,1:types(tt).nNodes),'rows') ;
        % Process slice configurations
            for cc = 1:max(sliceConfig)
            % Get corresponding elements
                isConf = sliceConfig==cc ;
                if ~any(isConf) ; continue ; end
            % Indices of nodes and edges of these elements 
                confIdx = nodeIdx(elmtIdx(isConf),[1:types(tt).nNodes nMaxNodesByElem+(1:types(tt).nEdges)]) ;
            % Fill new elements
                for pp = {'IN','ON','OUT'}
                    part = pp{1} ;
                    % Frome element indices to mesh indices
                        globIdx = S(cc).(part).dataAtIndices(confIdx') ;
                    % Type indices
                        ti = repmat(S(cc).(part).TypeIdx,[size(globIdx,3) 1]) ;
                    % Node indices
                        ni = reshape(permute(globIdx,[1 3 2]),[],size(globIdx,2)) ;
                    % Element table
                        newElems = pkg.geometry.mesh.elements.ElementTable('Types',S(cc).(part).Types,'Indices',[ti ni]) ;
                    sliceElems.(part) = [sliceElems.(part) newElems] ;
                end
            end
        end
    % Add non-sliced elements to the list
        sliceElems.IN = [sliceElems.IN this.Elems.subpart(fBool==-1)] ;
        sliceElems.OUT = [sliceElems.OUT this.Elems.subpart(fBool==1)] ;
    % Compute edge intersections with the levelset
        uEdg = unique(edgesOnLvlSt(:)) ;
        Xc = NaN(this.nEdges,this.nCoord) ;
        Xc(uEdg,:) = edgCrossLvlSet(this,fcn,uEdg,tol,X) ;
    % All Nodes
        X = [this.Nodes ; Xc] ;
    % Create the meshes
        mesh = [] ;
        for pp = {'IN','ON','OUT'}
            part = pp{1} ;
            mesh.(part) = pkg.geometry.mesh.Mesh('Nodes',X,'Elems',sliceElems.(part)) ;
            mesh.(part).cullUnused ;
        end
    end
    
    function mesh = split(this,fcn)
    % Split the mesh at locations where the sign of fcn(X) changes
    end
end


%% MESH CLEANUP
methods

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
        keep = ~remove ;
        if all(keep) ; return ; end
    % Remove nodes
        mesh.Nodes = mesh.Nodes(keep,:) ;
    % Change the element indices
        newNodeIdx = double(keep) ;
        newNodeIdx(keep) = 1:sum(keep) ;
        mesh.Elems = mesh.Elems.changeNodeIdx(newNodeIdx(:)) ;
    end

    function mesh = removeElems(this,remove)
    % Remove elems in the list and clean the mesh
    % nod can be logical or indexes
        if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
        if ~islogical(remove) % convert to logical
            remove = full(sparse(remove(:),ones(numel(remove),1),true,mesh.nElems,1)) ; 
        end
        if ~any(remove) ; return ; end
        mesh.Elems = mesh.Elems.subpart(~remove) ;
        mesh.cullUnused ;
    end

    function iNodes = invalidNodes(this)
    % Detect invalid nodes (Inf or NaN)
        iNodes = any( isinf(this.Nodes) | isnan(this.Nodes) ,2) ;
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
        [~,new,old] = uniquetol(mesh.Nodes,tol,'ByRows',true,'DataScale', 1) ;
    % Take the barycenter for duplicates groups
        meanMat = sparse(old,1:numel(old),1,numel(new),numel(old)) ;
        meanMat = sparse(1:numel(new),1:numel(new),1./sum(meanMat,2))*meanMat ;
    % Return moved node if needed
        if nargout>1 ; nodesMoved = mesh.Nodes(any(mod(meanMat,1),1),:) ; end
    % Change the node list
        mesh.Nodes = meanMat*mesh.Nodes ;
    % New element list
        mesh.Elems = mesh.Elems.changeNodeIdx(old) ;
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


%% EXPERIMENTAL
methods
    function mesh = fillHoles(this,method)
    % Return a mesh that corresponds to the filling of holes in the mesh
    % works in 2D only for now
        if nargin<2 ; method = 'distmesh' ; end
        mesh = copy(this) ; % always take a copy
        X = this.Nodes ;
    % Extract boundary curves
        [~,bndCrv] = this.boundaryCurves ;
    % Create new meshes from the inner contours and merge with the current mesh
        for bb = 2:numel(bndCrv) % the outer boundary curve is often the first; we skip it
            switch method
                case 'distmesh'
                    x = X(bndCrv{bb},:) ; 
                    x = unique(x,'rows','stable') ;
                    lvlst = pkg.geometry.levelset.Polygon(x) ;
                    dens = pkg.geometry.density.Polygon(x) ;
                    h0 = min(dens.Le) ;
                    fh = @dens.evalAt ;
                    fillMesh = pkg.geometry.mesh.distMesh(lvlst...
                                    ,'h0',h0 ...
                                    ,'fh',fh ...
                                    ,'t_dmax',0 ...
                                    ,'p_dmax',1e-6 ...
                                    ,'showMesh',false ...
                                ) ;
                otherwise
                    fillMesh = pkg.geometry.mesh.fromBoundary(X(bndCrv{bb},:)) ;
            end
            mesh = merge(mesh,fillMesh) ;
        end
    end
end

%% ELEMENT CONVERSION
methods
    function mesh = simplex(this)
    % Convert to a simplex mesh
        mesh = copy(this) ;
        mesh.Elems = mesh.Elems.simplex ;
    end
    
    function mesh = quadhex(this)
    % Convert as much elements as possible to quad/hex elements
    % First convert to simplex (tri/tetra elems)
        mesh = simplex(this) ;
    % Triangles to quads
        % Find triangles
            isTri = isa(mesh.Elems.Types,'pkg.geometry.mesh.elements.base.Triangle') ;
            isTri = ismember(mesh.Elems.TypeIdx,find(isTri)) ;
            tri = find(isTri) ;
        % Attached interior edges
            ele2edg = mesh.elem2edge ; 
            intE = sum(ele2edg(:,tri),2)==2 ;
            nE = full(sum(intE)) ;
        % Attached triangles
            [ee,tt] = find(ele2edg(intE,tri)) ;
            [~,is] = sort(ee) ;
            tt = reshape(tt(is),2,[])' ; % [nE 2]
        % [Edg setdiff(Tri,Edg)] Nodes
            nn = mesh.Elems.NodeIdx(tri(tt),1:3) ; % [nE*2 3]
            nn = [mesh.Edges.NodeIdx(intE,1:2) reshape(nn,nE,6)] ; % [nE 8]
        % Quads : keep only tri nodes that are not edg nodes
            nn = reshape(real(unique(double(nn)'+(1:nE)*1i,'stable')),[4 nE])' ;
            quads = nn(:,[1 3 2 4]) ; % [nE 4]
        % Quad angles: scalar product of edge vectors
            xq = reshape(mesh.Nodes(quads,:),nE,4,mesh.nCoord) ;
            vec = xq(:,[2:end 1],:)-xq ; % [nE 4 mesh.nCoord]
            vec = vec./sqrt(sum(vec.^2,3)) ; % [nE 4 mesh.nCoord]
            scal = sum(vec.*vec(:,[2:end 1],:),3) ; % % [nE 4] should be in [0 1]
        % Quad quality associated to the maximum angle deviation from 90°
            Q = 1-max(abs(scal),[],2) ; % [nE 1]
        % Sort edges by quad quality
            [~,ee] = sort(Q,'descend') ;
            tt = tt(ee,:) ;
            quads = quads(ee,:) ;
        % Cull edges that are attached to triangles already used
            ee = [] ; ttt = unique(tt(:)) ;
            for eee = 1:nE
            % all the triangles attached to the edge are still available ?
                if ~all(ismember(tt(eee,:),ttt)) ; continue ; end
            % Add the edge
                ee = [ee eee] ;
            % Remove the edge's triangles from the list
                ttt = setdiff(ttt,tt(eee,:)) ;
                if isempty(ttt) ; break ; end
            end
        % Create the quad elements
            qElmt = pkg.geometry.mesh.elements.base.Quadrangle ;
            qIdx = padarray(quads(ee,:),[0 1],1,'pre') ;
            quads = pkg.geometry.mesh.elements.ElementTable('Types',qElmt,'Indices',qIdx) ;
        % Retrieve triangles that ave not been joined
            tri = tri(setdiff(1:numel(tri),tt(ee,:))) ;
    % Join all elements
        mesh.Elems = [mesh.Elems.subpart(tri) quads] ;
    % Sort elements
        mesh.sortElems ;
    end
end

%% VIZUALIZATION
methods
    function h = plot(this,varargin)
    % Simple plot of the mesh
        h = pkg.geometry.mesh.MeshPlot(varargin{:}) ;
        h.Parent = gca ; 
        h.Mesh = this ;
    end
end


%% IMPORT & EXPORT FUNCTIONS
methods
    function stlwrite(this,filename)
    % Write an STL file from the mesh
        if nargin<2
            [file,path] = uiputfile('*.stl','Export the mesh to STL','mesh') ;
            if path==0 ; return ; end
            filename = [path file] ;
        end
    % Take only the outer surface
        faces = this.Faces.subpart(this.outerFaces) ;
    % Convert to triangle mesh
        faces = faces.simplex ;
    % Cull unused nodes
        nodeIdx = faces.NodeIdx ;
        [uIdx,~,ic] = unique([nodeIdx(:) ; 0],'sorted') ; % force the zero to be there
        uIdx = uIdx(2:end) ; % remove the zero index
        nodeIdx = reshape(ic(1:end-1)-1,size(nodeIdx)) ;
    % Create a MATLAB-builtin triangulation
        tri = triangulation(double(nodeIdx),this.Nodes(uIdx,:)) ;
    % Write the file
        stlwrite(tri,filename) ;
    end
end




end

