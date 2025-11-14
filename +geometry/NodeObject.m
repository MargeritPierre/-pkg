classdef (Abstract) NodeObject < matlab.mixin.Copyable
%NODEOBJECT (Abstract) Interface for a geometry defined by nodes
% Implements several useful methods associated to a nodal definition
    
%% ABSTRACT PROPERTIES & METHODS
% That have to be instanciated in inherited subclasses
methods (Abstract)
    normals = nodeNormals(this)
    % Return the normal to the geometry evaluated on nodes
    % normals [nNodes nCoord]
    
    obj = replicate(this,x,fill,close)
    % Replication function given a 3D array of coordinates x representing a
    % trajectory x [nNodes nCoord nPtsOnTrajectory]
    % used in extrusion, sweep, revolution, etc..
end


%% PROPERTIES
properties (SetObservable)
% The list of nodes [nNodes nCoord]
    Nodes
end
    

%% CONSTRUCTOR / Destructor
methods
    function this = NodeObject(varargin)
    % Constructor
    end

    function delete(this)
    % Destructor
    end
end


%% COUNT FUNCTIONS
methods
    function N = nNodes(this) ; [N,~] = cellfun(@size,{this.Nodes}) ; end
end
properties (Dependent) 
% The number of coordinates of each node
    nCoord
end


%% GEOMETRICAL DESCRIPTION
methods
% NUMBER OF COORDINATES
% Allow to easily change the mesh space dimension (any integer>0)
    function nCoord = get.nCoord(this) 
    % Get the numer of coordinates
        nCoord = size(this.Nodes,2) ; 
    end
    function set.nCoord(this,nCoord)
    % Set the number of coordinates
        this.Nodes = [this.Nodes(:,1:min(this.nCoord,nCoord)) zeros(this.nNodes,nCoord-this.nCoord)] ;
    end
    
% APPROXIMATE SIZE (Bounding Box) [2 nCoord] (min ; max)
    function bbox = boundingBox(this)
        bbox = [min(this.Nodes,[],1) ; max(this.Nodes,[],1)] ;
    end
    
% PLANARITY
    function [ispl,normal,origin,frame] = isPlanar(this,tol)
    % Check for the object planarity
    % Also return the plane normal, origin 
    % and full frame [ivect icoord]
        origin = mean(this.Nodes,1,'omitnan') ;
        % If the object is trivially plane
            if this.nCoord<3 % 1D or 2D coordinates
                ispl = true ; 
                normal = [0 0 1] ;
                frame = eye(3) ;
                return ;
            end
        % Otherwise, try to find a common plane by finding the
        % principal coordinates variances
            if nargin<2 ; tol = this.defaultTolerance ; end
            % Compute the coordinates covariance
                P = this.Nodes-origin ;
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
    
% CLOSEST NODE
    function [nn,dist] = closestNode(this,P,X)
    % Return the closest node indices
    % P = [nPts nCoord] or [nPts 1 nCoord]
    % P = [nNodes nCoord] or [nNodes 1 nCoord]
    % nn = [nPts 1]
        if nargin<3 ; X = this.Nodes ; end
        [nn,sqdist] = pkg.data.closestPoint(P,X) ;
        if nargout>1 ; dist = sqrt(sqdist) ; end
    end
end


%% GEOMETRICAL OPERATIONS
methods 
    function tol = defaultTolerance(this,x)
    % Return the default tolerance used in geometrical operations
    % Optional X [nNodes nCoord] can be used for custom node coordinates
        if nargin<2 ; x = this.Nodes; end
        tol = norm(range(x,1))*1e-9 ;
    end

    function obj = plus(A,B)
    % Addition of two objects
        obj = copy(A) ; % Always take a copy
        % New mesh coordinates values
            nNodes = A.nNodes + B.nNodes ;
            x = NaN(nNodes,max(A.nCoord,B.nCoord)) ;
            x(1:A.nNodes,1:A.nCoord) = A.Nodes ;
            x(A.nNodes+1:end,1:B.nCoord) = B.Nodes ;
        % Return/Modify the object
            obj.Nodes = x ; 
    end
    
    function [obj,nodA,nodB] = merge(A,B,tol,extBox)
    % Merge two meshes
    % nodA/nodB contains the indices of the output nodes that correspond to
    % each node in A/B
    % obj is the same type than A
        if nargin<3 ; tol = A.defaultTolerance ; end
        if nargin<4 ; extBox = 1.2 ; end % allow to extend the merging bbox
    % Find the closest nodes
        % Reduce to bounding boxes
            bboxA = [min(A.Nodes,[],1) ; max(A.Nodes,[],1)] + extBox*[-1 ; 1]*tol ;
            bboxB = [min(B.Nodes,[],1) ; max(B.Nodes,[],1)] + extBox*[-1 ; 1]*tol ;
            AinB = find(all(A.Nodes>=bboxB(1,:) & A.Nodes<=bboxB(2,:),2)) ;
            BinA = find(all(B.Nodes>=bboxA(1,:) & B.Nodes<=bboxA(2,:),2)) ;
        % Compute the distances
            sqDist = sum((permute(B.Nodes(BinA,:),[1 3 2])-permute(A.Nodes(AinB,:),[3 1 2])).^2,3) ;
            [nodB,nodA] = find(sqDist<=tol.^2) ; % all pairs of close nodes
            sqDist = sqDist(sub2ind(size(sqDist),nodB,nodA)) ;
        % Convert to global numerotation (including out-of-bbox)
            nodA = AinB(nodA) ;
            nodB = BinA(nodB) ;
    % Do not merge multiple nodes with one node (!!!)
        % Sort by increasing distance
            [~,indSort] = sort(sqDist,'ascend') ;
            nodB = nodB(indSort) ; 
            nodA = nodA(indSort) ;
        % Take only the closest node that is not already taken
            [nodA,nn] = unique(nodA,'stable') ;
            nodB = nodB(nn) ;
            [nodB,nn] = unique(nodB,'stable') ;
            nodA = nodA(nn) ;
    % Merge Nodes
        newNodes = mean(cat(3,A.Nodes(nodA,:),B.Nodes(nodB,:)),3) ;
        % Complete list
            x = [A.Nodes ; B.Nodes] ;
            newNodeIdx = 1:size(x,1) ;
        % Replace merged nodes
            x(nodA,:) = newNodes ;
            x(A.nNodes + nodB,:) = newNodes ;
        % Cull duplicates
            newNodeIdx(A.nNodes + nodB) = nodA ;
            newNodeIdx(setdiff(1:numel(newNodeIdx),A.nNodes+nodB)) = 1:numel(newNodeIdx)-numel(nodB) ;
            x(A.nNodes + nodB,:) = [] ;
    % Keep the lists of replacement nodes
        nodA = newNodeIdx(1:A.nNodes) ;
        nodB = newNodeIdx(A.nNodes+(1:B.nNodes)) ;
    % Return the new object
        obj = copy(A) ; obj.Nodes = x ;
    end

    function obj = applyTransform(this,F,v)
    % Apply a geometrical transformation to the object
    % F is the transfo: [nCoord nCoord N] array,
    % v is the translation vector; has to be of size [N nCoord]
    % with N>=1 to replicate the object on different transfos.
        if nargout==0 ; obj = this ; else ; obj = copy(this) ; end
        if nargin<3 ; v = zeros(1,size(F,2)) ; end
        % Force formatting
            F = F(:,:,:) ; v = permute(v(:,:),[3 2 1]) ; % [1 nCoord N]
            F = F + v*0 ; v = v + F(1,:,:)*0 ;
            [nC,~,N] = size(F) ;
        % Fix incompatible dimensions
            dC = obj.nCoord-nC ;
            if dC>0
                O = zeros(nC,dC,N) ;
                F = [F O ; permute(O,[2 1 3]) repmat(eye(dC),[1 1 N])] ;
                v = [v O(1,:,:)] ;
            elseif dC<0
                obj.nCoord = nC ;
            end
        % Apply the transform to the coordinates
            x = pkg.math.innerprod(this.Nodes,F) + v ;
        % Replicate the object (without filling)
            obj.replicate(x,false) ;
    end

    function obj = move(this,T)
    % Move the object by a translation vector T
    % T has to be of size [N nCoord] with N>=1 to replicate the mesh
        if nargout==0 ; obj = this ; else ; obj = copy(this) ; end
        % Format the translation vector
            T = squeeze(T) ;
        % Apply transform
            obj.applyTransform(eye(size(T,2)),T) ;
    end

    function obj = rotate(this,THETA,CEN,AX)
    % Rotate the object by 
    % the angle THETA: [N 1]  
    % around the point P: [1 nCoord] (default [0 0 0])
    % with respect to the axis AX: [1 3]
    % with N>=1 to replicate the object
        if nargout==0 ; obj = this ; else ; obj = copy(this) ; end
        if nargin<3 ; CEN = zeros(1,obj.nCoord) ; end
        if nargin<4 ; AX = [0 0 1] ; end
        % Move the mesh
            obj.move(-CEN) ;
        % Compute the transformation matrix
            F = pkg.math.rotmat(-THETA,AX) ;
            if this.nCoord<3 && all(all(AX==[0 0 1],2))
                F = F(1:2,1:2,:) ;
            end
        % Rotate
            obj.applyTransform(F) ;
        % Re-move the mesh
            obj.move(CEN) ;
    end

    function obj = scale(this,SC,CEN)
    % Scale the object by 
    % the scale factors: [N nCoord]  
    % w.r.t the center point P: [N nCoord] (default [0 0 0])
    % with N>=1 to replicate the object
        if nargout==0 ; obj = this ; else ; obj = copy(this) ; end
        if nargin<3 ; CEN = zeros(1,obj.nCoord) ; end
    % Match argument size
        SC = permute(SC,[3 2 1]) ; CEN = permute(CEN,[3 2 1]) ;
        [x,CEN] = pkg.data.matchsize(2,obj.Nodes,CEN,'filler',0) ;
        if size(SC,2)==1 % isotropic scaling
            SC = repmat(SC,1,size(x,2)) ;
        else
            SC(:,end+1:size(x,2),:) = 1 ; % unit scale by default
        end
    % Apply the transform
        x = (x-CEN).*SC + CEN ;
    % Replicate the object (without filling)
        obj.replicate(x,false) ;
    end
    
    function obj = offset(this,DIST,fill)
    % Offset the object
    % DIST contains the offset distances (can be negative)
    % fill: (bool) set true when you want to create..
    % ..a surface/volume from curves/surfaces
    % fill: bool. false: just offset the object; true: fill the output
        if nargout==0 ; obj = this ; else ; obj = copy(this) ; end
    % Process inputs
        if numel(DIST)==0 ; return ; end % return the same object
        if nargin<3 ; fill = false ; end % do not fill by default
        if fill && numel(DIST)==1 ; DIST = [0 ; DIST(:)] ; end % force the volume to be valid
        DIST = sort(DIST(:),'descend') ; 
        nOfst = numel(DIST) ;
    % Offset nodes
        nodeNormals = this.nodeNormals ;
        x = this.Nodes + reshape(DIST,[1 1 nOfst]).*nodeNormals(:,1:this.nCoord) ;
    % Replicate elements with the given coordinates array
        obj.replicate(x,fill) ;
    end
    
    function obj = sweep(this,CRV,mode)
    % Create a sweep mesh geometry from an object and a curve
    % CRV: [nPts nCoord]
    % mode: 'angle' or 'parallel' (default)
    %   - 'angle': the replicated objects remain at the same angle w.r.t the curve
    %   - 'parallel': all nodes are swept by the same curve (parallel edges)
        if nargout==0 ; obj = this ; else ; obj = copy(this) ; end
    % Process inputs
        if nargin<3 ; mode = 'parallel' ; end
    % Initialize
        x = this.Nodes ;
    % Match Curve & nodes dimensions
        nC = max(size(x,2),size(CRV,2)) ; % target nCoord
        x = [x zeros(size(x,1),nC-size(x,2))] ;
        CRV = [CRV zeros(size(CRV,1),nC-size(CRV,2))] ; % [nPts nC]
    % Sweep nodes
        switch mode
            case 'parallel' % x = x + u0
                u0 = permute(CRV-CRV(1,:),[3 2 1]) ; % [1 nC nPts] ;
                x = x + u0 ; % [nNodes nC nPts]
            case 'angle' % x = Rotate(x-crv(1),angle(tang_1,tang_i)) + crv
                % Tangent vectors [nPts nC]
                    tang = diff(CRV,1,1) ; % segments
                    tang = tang./sqrt(sum(tang.^2,2)) ; % normalize
                    tang = [3*tang(1,:)-tang(2,:) ; tang(1:end-1,:)+tang(2:end,:) ; 3*tang(end,:)-tang(end-1,:)] ; % node tangents
                    tang = tang./sqrt(sum(tang.^2,2)) ; % normalize
                % Rotation vector R = vectprod(t_1,t_n)
                    % R = pkg.math.vectprod(tang(1,:),tang,2) ; % [nPts nC]
                    R = pkg.math.vectprod(tang([1,1:end-1],:),tang,2) ; % [nPts nC]
                % Rotation matrices
                    F = pkg.math.rotmat(-R) ; % [3 3 nPts]
                    for pp = 2:size(F,3)
                        F(:,:,pp) = F(:,:,pp-1)*F(:,:,pp) ;
                    end
                    if nC<3 ; F = F(1:nC,1:nC,:) ; end
                % Apply
                    x = pkg.math.innerprod(x-CRV(1,:),F) + permute(CRV,[3 2 1]) ;
            otherwise
                error('Unknown sweeping mode')
        end
    % Replicate elements given the new coordinates array
        obj.replicate(x,true) ;
    end
    
    function obj = extrude(this,VEC,DIV)
    % Extrude the object along the vector VEC
    % (opt.) use DIV (int. scalar) to subdivide the extrusion
    % VEC: [1 nCoord]
        if nargin<3 ; DIV = 1 ; end
        if nargout==0 ; obj = this ; else ; obj = copy(this) ; end
        VEC = VEC(:)' ;
    % Match coordinates
        obj.nCoord = max(obj.nCoord,size(VEC,2)) ;
        VEC = padarray(VEC,[0 obj.nCoord-size(VEC,2)],0,'post') ;
    % Sweeping curve
        CRV = interp1(VEC.*[0;1],linspace(1,2,DIV+1)) ;
        x = obj.Nodes + permute(CRV,[3 2 1]) ;
        obj.replicate(x,true) ; 
    end
    
    function obj = revolution(this,AX,PT,ANG,DIV,loop) 
    % Revolve an object around AX centered on PT with angle ranges ANG
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
        if nargin<6 ; loop = abs(ANG(end)-ANG(1)-2*pi)<100*eps ; end
        if nargout==0 ; obj = this ; else ; obj = copy(this) ; end
    % Rotation matrices
        R = pkg.math.rotmat(-ANG(:),AX) ;
    % New coordinates
        x = pkg.math.innerprod(obj.Nodes-PT,R,2,1) + PT ;
    % Apply
        obj.replicate(x,true,loop) ;
    end   
    
end

end

