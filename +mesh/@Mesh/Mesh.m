classdef Mesh < matlab.mixin.Copyable
% The general mesh class

    
%% CONSTRUCTOR / DESTRUCTOR
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
   

%% CONNECTIVITY LISTS
    methods
        
        function [edges,ie] = listEdges(this,list)
        % Return the unique edges associated to a node list and the list of
        % corresponding edges indices
            if nargin<2 ; list = this.Elems ; end
            % Number of nodes in each element
                Nn = sum(~isnan(list),2) ;
            % Edge nodes
                n1 = list ;
                n2 = circshift(list,1,2) ;
            % Replace NaNs by the correct index
                list(:,1) = n2(sub2ind(size(n2),1:size(list,1),Nn(:)')) ;
            % All edges
                edges = [n1(:) n2(:)] ;
            % NaN edges
                nans = any(isnan(edges),2) ;
                edges(nans,:) = [] ;
            % Unique edges
                edges = sort(edges,2) ;
                [edges,~,iu] = unique(edges,'rows') ;
            % Reshape edges indices
                ie = NaN(size(list)) ;
                ie(~nans) = iu ;
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
        
        function [crv,asCell] = boundaryCurves(this)
        % Return an array of node indices representing the sorted mesh boundaries
        % crv contains the nodes indices with the different curves
        % separated by NaNs
            if this.nElems==0 ; crv = {} ; return ; end
            % List of boundary edges
                bndEdg = boundaryEdges(this) ;
                nn = this.Edges.NodeIdx(bndEdg,:) ;
            % Scan the edges to build the curves
                crv{1} = nn(1,:) ;
                nn(1,:) = [] ;
                while ~isempty(nn)
                    % Find the next point
                        [nextEdg,indNode] = find(nn==crv{end}(end)) ;
                    % If no next point found, create a new curve
                        if isempty(nextEdg)
                            crv{end+1} = nn(1,:) ;
                            nn(1,:) = [] ;
                    % Otherwise, add it
                        else
                            crv{end} = [crv{end} nn(nextEdg(1),3-indNode(1))] ;
                            nn(nextEdg(1),:) = [] ;
                        end
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
        
%         function P = list2Points(this,list,P)
%         % Convert a list of node indices to the associated coordinates
%         % size(P) = [size(list) nCoord]
%             if nargin<3 ; P = this.Nodes ; end
%             nans = isnan(list) ;
%             list(nans) = 1 ;
%             P = P(list(:),:) ;
%             if islogical(P) ; P(nans,:) = false ;
%             else ; P(nans,:) = NaN ; end
%             P = reshape(P,[size(list) size(P,2)]) ;
%         end
        
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
        end
    end
    
%% MESH OPERATIONS (MOTION, BOOLEAN)
    methods 
        
        function mesh = plus(this,mesh2)
        % Addition of two meshes (NEED CLEANUP TO WELD THE MESHES!)
            if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end
            % New mesh coordinates values
                nNodes = this.nNodes + mesh2.nNodes ;
                x = NaN(nNodes,max(this.nCoord,mesh2.nCoord)) ;
                x(1:this.nNodes,1:this.nCoord) = this.X.Values ;
                x(this.nNodes+1:end,1:mesh2.nCoord) = mesh2.X.Values ;
            % New element table
                elems = [this.Elems mesh2.Elems] ;
                elems.NodeIdx(this.nElems+1:end) = elems.NodeIdx(this.nElems+1:end) + this.nNodes ;
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
        
        function tol = defaultTolerance(this)
        % Return the default tolerance used in cleanup operations
            tol = norm(range(this.X.Values))*1e-6 ;
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
            iElems = any(this.Elems.NodeIdx>this.nNodes) ; % At least one node index too large
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
            mesh.Elems = mesh.list2Points(mesh.Elems,old);
        end
        
    end
    
    
%% SUBDIVISION, SMOOTHING
    methods % Functions defined in external files
        % Laplacian Smoothing
        mesh = LaplacianSmoothing(this,lmbda,iterations,tol)
        % Mesh Subdivision
        mesh = CatmullClark(this,iter)
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

