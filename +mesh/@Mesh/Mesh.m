classdef Mesh < matlab.mixin.Copyable
% The general mesh class

    
%% CONSTRUCTOR / DESTRUCTOR
    methods
        function this = Mesh(varargin)
        % Class Constructor
            if mod(nargin,2) ; error('wrong number of arguments') ; end
        % Initialize the node position
            this.X = pkg.mesh.MeshFunction('Mesh',this) ;
        % Process input arguments
            for arg = 1:2:nargin-1
                this.(varargin{arg}) = varargin{arg+1} ;
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
                    if ~ismatrix(elems) ; error('wrong format: must be [nElems maxNodesByElems]') ; end
                    elems = uint32(elems) ;
                    elems = pkg.mesh.elements.ElementTable('Indices',elems,'Types',this.Elems.Types) ;
            end
        % Clean the element list
            elems = clean(elems) ;
        % Set mesh features
            this.Faces = elems.uniqueFaces ;
            this.Edges = elems.uniqueEdges ;
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
        function N = nNodes(this) ; N = max(this.Elems.Indices(:)) ; end
        function N = nElems(this) ; N = size(this.Elems.Indices,1) ; end
        function N = nEdges(this) ; N = size(this.Edges.Indices,1) ; end
        function N = nFaces(this) ; N = size(this.Faces.Indices,1) ; end
        % Other ...
        function N = maxNodesByElem(this) ; N = size(this.Elems.Indices,2) ; end
        function N = nNodesInElem(this) ; N = this.Elems.nNodes ; end
    end
    
%% MESH SPACE DIMENSION
    properties (Dependent) ; nCoord ; end
    methods
        % Allow to easily change the mesh space dimension (any integer>0)
        function set.nCoord(this,nCoord)
            this.X = [this.X.Values(:,1:min(this.nCoord,nCoord)) zeros(this.nNodes,nCoord-this.nCoord)] ;
        end
        function nCoord = get.nCoord(this) ; nCoord = size(this.X) ; end
    end
   

% %% CONNECTIVITY LISTS
%     methods
%         function this = set.Elems(this,elems)
%         % When the element list change, perform some checks and
%             % Format check
%                 if ~ismatrix(elems) ; error('wrong format: must be [nElems maxNodesByElems]') ; end
%             % Type recasting : stricly positive integer or NaN
%                 elems = round(elems) ;
%                 elems(elems<1) = NaN ;
%             % Cull duplicates
%                 elems(isnan(elems)) = 0 ;
%                 [~,ie,~] = unique(sort(elems,2),'rows') ;
%                 elems = elems(ie,:) ;
%                 elems(elems==0) = NaN ;
%             % Set
%                 this.Elems = elems ;
%             % Compute other connectivity lists
%                 this = this.setFaces() ;
%                 this = this.setEdges() ;
%             % Add dummy nodes if needed
%                 if max(elems(:))>this.nNodes
%                     %warning('adding dummy nodes') ;
%                     %this.Nodes(this.nNodes+1:max(elems(:)),:) = NaN ;
%                 end
%         end
%         
%         function this = setFaces(this)
%         % Compute the faces from the element list
%             if this.nElems==0 ; this.Faces = [] ; return ; end
%             this.Faces = this.Elems ; % <TODO> How to differentiate faces and elements ?
%         end
%         
%         function this = setEdges(this)
%         % Compute the edges from the element list
%             if this.nElems==0 ; this.Edges = [] ; return ; end
%             % Unique edges
%                 edges = listEdges(this) ;
%             % Set
%                 this.Edges = edges ;
%         end
%         
%         function [edges,ie] = listEdges(this,list)
%         % Return the unique edges associated to a node list and the list of
%         % corresponding edges indices
%             if nargin<2 ; list = this.Elems ; end
%             % Number of nodes in each element
%                 Nn = sum(~isnan(list),2) ;
%             % Edge nodes
%                 n1 = list ;
%                 n2 = circshift(list,1,2) ;
%             % Replace NaNs by the correct index
%                 list(:,1) = n2(sub2ind(size(n2),1:size(list,1),Nn(:)')) ;
%             % All edges
%                 edges = [n1(:) n2(:)] ;
%             % NaN edges
%                 nans = any(isnan(edges),2) ;
%                 edges(nans,:) = [] ;
%             % Unique edges
%                 edges = sort(edges,2) ;
%                 [edges,~,iu] = unique(edges,'rows') ;
%             % Reshape edges indices
%                 ie = NaN(size(list)) ;
%                 ie(~nans) = iu ;
%         end
%             
%     end
%    
%     
% %% CONNECTIVITY MATRICES
%     methods
%         function M = list2nod(this,list,varargin)
%         % Return a connectivity list deduced from a list of node indices
%             % Default values format
%                 if nargin<3 
%                     values = 'logical' ; 
%                 else
%                     values = varargin{1} ;
%                 end
%             % Switch values format
%                 switch values
%                     case 'logical'
%                         vvv = true(size(list)) ;
%                     case 'indices'
%                         vvv = repmat(1:size(list,2),[size(list,1) 1]) ;
%                     case 'mean'
%                         vvv = true(size(list)) ;
%                     otherwise
%                         vvv = varargin{1}(:) + list(:)*0 ;
%                 end
%             % Build the matrix
%                 iii = list(:) ;
%                 jjj = repmat(1:size(list,1),[1 4]) ;
%                 valid = ~isnan(iii) ;
%                 M = sparse(iii(valid),jjj(valid),vvv(valid),max(list(:)),size(list,1)) ; 
%             % Eventually compute the mean
%                 if strcmp(values,'mean') ; M = M./sum(M,2) ; end
%         end
%     end
% 
%     methods
%         function M = elem2nod(this,varargin)
%         % Elements to Nodes [nNodes nElems]: nod(:) = M*elmt(:)
%             M = this.list2nod(this.Elems,varargin{:}) ;
%         end
%         function M = face2nod(this,varargin)
%         % Faces to Nodes [nNodes nFaces]: nod(:) = M*fa(:)
%             M = this.list2nod(this.Faces,varargin{:}) ;
%         end
%         function M = edge2nod(this,varargin)
%         % Faces to Nodes [nNodes nFaces]: nod(:) = M*fa(:)
%             M = this.list2nod(this.Edges,varargin{:}) ;
%         end
%         function M = elem2edge(this,varargin)
%         % Elements to Edges [nEdges nElems]: edg(:) = M*elmt(:)
%         % An edge is linked to an element only if they share two nodes
%             %M = this.edge2nod'*this.elem2nod==2 ;
%             [~,ie] = this.listEdges(this.Elems) ;
%             M = this.list2nod(ie,varargin{:}) ;
%         end
%         function M = face2edge(this,varargin)
%         % Faces to Edges [nFaces nElems]: edg(:) = M*face(:)
%         % An edge is linked to a face only if they share two nodes
%             %M = this.edge2nod'*this.face2nod==2 ;
%             [~,ie] = this.listEdges(this.Faces) ;
%             M = this.list2nod(ie,varargin{:}) ;
%         end
%     end
%     
%     
% %% BOUNDARIES
%     methods
%         
%         function bndEdg = boundaryEdges(this)
%         % Return an array of logical, true if the edge is on the boundary
%         % The edge is on the boundary if it is linked to only one element
%             bndEdg = sum(this.elem2edge,2)==1 ;
%         end
%         
%         function bndNod = boundaryNodes(this)
%         % Return an array of logical, true if the node is on the boundary
%         % The node is on the boundary if it is linked to a boundary edge
%             bndEdg = boundaryEdges(this) ;
%             bndNod = logical(this.edge2nod*bndEdg) ;
%         end
%         
%         function [crv,asCell] = boundaryCurves(this)
%         % Return an array of node indices representing the sorted mesh boundaries
%         % crv contains the nodes indices with the different curves
%         % separated by NaNs
%             if this.nElems==0 ; crv = {} ; return ; end
%             % List of boundary edges
%                 bndEdg = boundaryEdges(this) ;
%                 nn = this.Edges(bndEdg,:) ;
%             % Scan the edges to build the curves
%                 crv{1} = nn(1,:) ;
%                 nn(1,:) = [] ;
%                 while ~isempty(nn)
%                     % Find the next point
%                         [nextEdg,indNode] = find(nn==crv{end}(end)) ;
%                     % If no next point found, create a new curve
%                         if isempty(nextEdg)
%                             crv{end+1} = nn(1,:) ;
%                             nn(1,:) = [] ;
%                     % Otherwise, add it
%                         else
%                             crv{end} = [crv{end} nn(nextEdg(1),3-indNode(1))] ;
%                             nn(nextEdg(1),:) = [] ;
%                         end
%                 end
%             % Concatenate curves
%                 crv = reshape(crv,1,[]) ;
%                 asCell = crv ; % Backup individual curves
%                 crv(end+1,:) = {NaN} ;
%                 crv = cat(2,crv{:}) ;
%                 crv = crv(1:end-1) ;
%         end
%         
%     end
%     
%     
% %% GEOMETRY
%     methods
%         
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
%         
%         function c = circumCenters(this,list)
%         % Return the circumcenters associated to the given list
%         % By default, list is this.Elems
%             if nargin<2 ; list = this.Elems ; end
%             % Get the list of coordinates corresponding to the list
%                 P = this.list2Points(list) ;
%             % Compute the mean
%                 c = mean(P,2,'omitnan') ;
%                 c = reshape(c,size(list,1),[]) ;
%         end
%         
%         function [ispl,normal,origin,frame] = isPlanar(this,tol)
%         % Check for the mesh planarity
%         % Also return the plane normal, origin 
%         % and full frame [ivect icoord]
%             origin = mean(this.Nodes,1) ;
%             % If the mesh is trivially plane
%                 if this.nCoord<3 % 1D or 2D coordinates
%                     ispl = true ; 
%                     normal = [0 0 1] ;
%                     frame = eye(3) ;
%                     return ;
%                 end
%             % Otherwise, try to find a common plane by finding the
%             % principal coordinates variances
%                 if nargin<2 ; tol = 1e-9 ; end
%                 % Compute the coordinates covariance
%                     P = this.Nodes-origin ;
%                     [V,s] = eig(P'*P,'vector') ;
%                 % V is the frame basis, s the variances
%                     [s,ind] = sort(s,'descend') ;
%                     V = V(:,ind) ;
%                 % Test the planarity
%                     ispl = s(end)/s(1)<tol ;
%                 % Return objects
%                     frame = V' ;
%                     normal = frame(end,:) ;
%         end
%         
%         function [normals,frames] = faceNormals(this)
%         % Return an array containing the face normals
%         % NORMALS size: [nFaces 3]
%         % Also return the face frames if needed
%         % FRAMES size: [nFaces 3 3] ([iface ivect icoord])
%             % First check if the mesh is planar
%                 [ispl,normal,~,frame] = isPlanar(this) ;
%                 if ispl 
%                     normals = repmat(normal(:)',[this.nFaces 1]) ;
%                     frames = repmat(reshape(frame,[1 3 3]),[this.nFaces 1 1]) ;
%                     return
%                 end
%             % Otherwise fit a plane on each face
%                 normals = zeros(this.nFaces,3) ;
%                 frames = zeros(this.nFaces,3,3) ;
%                 C = this.circumCenters(this.Faces) ;
%                 for ff = 1:this.nFaces
%                     ii = this.Faces(ff) ;
%                     ii = ii(~isnan(ii)) ;
%                     P = this.Nodes(ii,:) - C(ff) ;
%                     % Compute the coordinates covariance
%                         [V,s] = eig(P'*P,'vector') ;
%                     % V is the frame basis, s the variances
%                         [~,ind] = sort(s,'asccend') ;
%                         V = V(:,ind) ;
%                     % Return the objects
%                         frames(ff,:,:) = V' ;
%                         normals(ff,:) = frames(ff,end,:) ;
%                 end
%         end
%         
%         function this = sortElems(this,dir)
%         % Sort the element nodes of a planar mesh in clockwise 
%         % or counter-clockwise(==trigo) direction
%             if ~this.isPlanar ; error('Cannot sort nodes of a non-planar mesh') ; end
%             if nargin<2 ; dir = 'trigo' ; end
%             % Barycentered coordinates
%                 X = this.list2Points(this.Elems,this.Nodes) ;
%                 X = X-mean(X,2,'omitnan') ;
%             % Angle
%                 aaa = angle(X(:,:,1)+1i*X(:,:,2)) ;
%             % Sort indices
%                 [~,indSort] = sort(aaa,2) ;
%                 switch dir
%                     case 'trigo'
%                     case 'counter-clockwise'
%                     case 'clockwise'
%                         indSort = flip(indSort,2) ;
%                 end
%             % To linear indices
%                 indSort = sub2ind(size(this.Elems),repmat((1:this.nElems)',[1 size(this.Elems,2)]),indSort) ;
%             % New element list
%                 this.Elems = this.Elems(indSort) ;
%         end
%     end
%     
% %% MESH OPERATIONS (MOTION, BOOLEAN)
%     methods 
%         
%         function mesh = plus(mesh1,mesh2)
%         % Addition of two meshes (NEED CLEANUP TO WELD THE MESHES!)
%             % New node list
%                 nCoord = max(mesh1.nCoord,mesh2.nCoord) ;
%                 nNodes = mesh1.nNodes + mesh2.nNodes ;
%                 nodes = NaN(nNodes,nCoord) ;
%                 nodes(1:mesh1.nNodes,1:mesh1.nCoord) = mesh1.Nodes ;
%                 nodes(mesh1.nNodes+1:end,1:mesh2.nCoord) = mesh2.Nodes ;
%             % New element list
%                 nElems = mesh1.nElems + mesh2.nElems ;
%                 maxNodesByElem = max(mesh1.maxNodesByElem,mesh2.maxNodesByElem) ;
%                 elems = NaN(nElems,maxNodesByElem) ;
%                 elems(1:mesh1.nElems,1:mesh1.maxNodesByElem) = mesh1.Elems ;
%                 elems(mesh1.nElems+1:end,1:mesh2.maxNodesByElem) = mesh2.Elems + mesh1.nNodes ;
%             % Return the new mesh
%                 mesh = pkg.mesh.Mesh('Nodes',nodes,'Elems',elems) ;
%         end
%         
%         function this = replicate(this,N)
%         % Reproduce the mesh N times
%             if N<2 ; return ; end
%             this.Nodes = repmat(this.Nodes,[N 1]) ;
%             this.Elems = repmat(this.Elems,[N 1]) + kron((0:N-1)',ones(this.nElems,1)) ;
%         end
%         
%         function this = applyTransform(this,F,v)
%         % Apply a geometrical transformation to the mesh
%         % F is the transfo: [nCoord nCoord N] array,
%         % v is the translation vector; has to be of size [N nCoord]
%         % with N>=1 to replicate the mesh on different transfos.
%             if nargin<3 ; v = zeros(1,size(F,2)) ; end
%             % Force formatting
%                 F = F(:,:,:) ; v = permute(v(:,:),[3 2 1]) ; % [1 nCoord N]
%                 F = F + v*0 ; v = v + F(1,:,:)*0 ;
%                 [nC,~,N] = size(F) ;
%             % Fix incompatible dimensions
%                 dC = this.nCoord-nC ;
%                 if dC>0
%                     O = zeros(nC,dC,N) ;
%                     F = [F O ; permute(O,[2 1 3]) O(1:dC,1:dC,:)] ;
%                     v = [v O(1,:,:)] ;
%                 elseif dC<0
%                     this.nCoord = nC ;
%                 end
%             % New element list if needed
%                 this.Elems = repmat(this.Elems,[N 1]) + this.nNodes*kron((0:N-1)',ones(this.nElems,1)) ;
%             % Add a fourth dimension to stack replicated meshes
%                 F = permute(F,[4 1 2 3]) ; % [1 nCoord nCoord N]
%                 v = permute(v,[4 1 2 3]) ; % [1 1 nCoord N]
%             % Apply Transform
%                 nodes = sum(this.Nodes.*F,2) + v ; % [nNodes 1 nCoord N]
%             % Reshape
%                 nodes = permute(nodes,[1 4 3 2]) ; % [nNodes N nCoord]
%                 this.Nodes = reshape(nodes,[this.nNodes*N this.nCoord]) ; % [nNodes*N nCoord]
%         end
%         
%         function this = move(this,T)
%         % Move the mesh by a translation vector T
%         % T has to be of size [N nCoord] with N>=1 to replicate the mesh
%             % Format the translation vector
%                 T = squeeze(T) ;
%             % Apply transform
%                 this = this.applyTransform(eye(size(T,2)),T) ;
%         end
%         
%         function this = rotate(this,THETA,CEN,AX)
%         % Rotate the mesh by
%         % the angle THETA: [N 1]  
%         % around the point P: [1 nCoord] (default [0 0 0])
%         % with respect to the axis AX: [1 3] (default [0 0 1], not implemented)
%         % with N>=1 to replicate the mesh
%             if nargin<3 ; CEN = zeros(1,this.nCoord) ; end
%             if nargin<4 ; AX = [0 0 1] ; end
%             % Move the mesh
%                 this = this.move(-CEN) ;
%             % Compute the transformation matrix
%                 THETA = reshape(THETA(:),1,1,[]) ;
%                 F = [cos(THETA) sin(THETA) ; -sin(THETA) cos(THETA)] ;
%                 if this.nCoord>2 ; F(end+1,end+1,:) = 1 ; end
%             % Rotate
%                 this = this.applyTransform(F) ;
%             % Re-move the mesh
%                 this = this.move(CEN) ;
%         end
%         
%     end
%     
% %% MESH CLEANUP
%     methods
%         
%         function tol = defaultTolerance(this)
%         % Return the default tolerance used in cleanup operations
%             tol = norm(range(this.Nodes))*1e-6 ;
%         end
%         
%         function this = clean(this,tol)
%         % Apply the following cleanup functions to the mesh
%             if nargin<2 ; tol = this.defaultTolerance ; end
%             this = this.cullInvalid ;
%             this = this.cullUnused ;
%             this = this.cullDuplicates(tol) ;
%         end
%         
%         function this = removeNodes(this,remove)
%         % Remove nodes in the list and clean the mesh
%         % nod can be logical or indexes
%             if ~islogical(remove) % convert to logical
%                 remove = full(sparse(remove(:),ones(numel(remove),1),true,this.nNodes,1)) ; 
%             end
%             if ~any(remove) ; return ; end
%             this.Nodes = this.Nodes(~remove,:) ;
%             this.Elems = this.list2Points(this.Elems,cumsum(~remove(:))) ;
%         end
%         
%         function this = removeElems(this,remove)
%         % Remove elems in the list and clean the mesh
%         % nod can be logical or indexes
%             if ~islogical(remove) % convert to logical
%                 remove = full(sparse(remove(:),ones(numel(remove),1),true,this.Elems,1)) ; 
%             end
%             if ~any(remove) ; return ; end
%             this.Elems = this.Elems(~remove,:) ;
%             this = this.cullUnused ;
%         end
%         
%         function iNodes = invalidNodes(this)
%         % Detect invalid nodes (Inf or NaN)
%             iNodes = any( isinf(this.Nodes) | isnan(this.Nodes) ,2) ;
%         end
%         
%         function iElems = invalidElems(this)
%         % Detect invalid elements (some nodes are invalid)
%             iElems = any(this.Elems>this.nNodes,2) ; % At least on node index too large
%             iElems = iElems | all(isnan(this.Elems),2) ; % All NaNs indices
%             % Points to invalid nodes
%             iNodes = invalidNodes(this) ;
%             iElems(~iElems) = any(this.list2Points(this.Elems(~iElems,:),iNodes(:)),2) ;
%         end
%         
%         function this = cullInvalid(this)
%         % Cull invalid nodes and elements
%             this = this.removeElems(this.invalidElems) ;   
%             this = this.removeNodes(this.invalidNodes) ;         
%         end
%         
%         function uNodes = unusedNodes(this)
%         % Return true for nodes not part of element list
%             uNodes = ~ismember(1:this.nNodes,this.Elems(:)) ;
%         end
%         
%         function this = cullUnused(this)
%         % Cull unused nodes
%             this.removeNodes(this.unusedNodes) ;
%         end
%         
%         function [this,nodesMoved] = cullDuplicates(this,tol)
%         % Cull duplicate nodes in the mesh and corresponding elements
%             if nargin<2 ; tol = this.defaultTolerance ; end
%         % Unique node indices list
%             [~,new,old] = uniquetol(this.Nodes,tol,'ByRows',true,'DataScale', 1) ;
%         % Take the barycenter for duplicates groups
%             meanMat = sparse(old,1:numel(old),1,numel(new),numel(old)) ;
%             meanMat = sparse(1:numel(new),1:numel(new),1./sum(meanMat,2))*meanMat ;
%         % Return moved node if needed
%             if nargout>1 ; nodesMoved = this.Nodes(any(mod(meanMat,1),1),:) ; end
%         % Change the node list
%             this.Nodes = meanMat*this.Nodes ;
%         % New element list
%             this.Elems = this.list2Points(this.Elems,old);
%         end
%         
%     end
%     
%     
% %% SUBDIVISION, SMOOTHING
%     methods % Functions defined in external files
%         % Laplacian Smoothing
%         this = LaplacianSmoothing(this,lmbda,iterations,tol)
%         % Mesh Subdivision
%         this = CatmullClark(this,iter)
%     end
%     
%     
% %% VIZUALIZATION
%     methods
%         function h = plot(this)
%         % Simple plot of the mesh
%             h = pkg.mesh.MeshPlot() ;
%             h.Parent = gca ; 
%             h.Mesh = this ;
%         end
%     end




end

