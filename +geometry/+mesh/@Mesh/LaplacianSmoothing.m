function mesh = LaplacianSmoothing(this,lmbda,iterations,tol)
%Apply Laplacian smoothing to a pkg.geometry.mesh.Mesh
% lmbda is a ratio btw 0 and 1; can be [lin lout], to give a different
% value for interior and boundary nodes
if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end

% Format inputs
    % Lambda
        if nargin<2 ; lmbda = [1 0]*0.5 ; end
        lmbda = lmbda(:)'.*[1 1] ;
    % Iterations
        if nargin<3 ; iterations = Inf ; end
        if nargin<4 ; tol = 1e-9 ; end
    % Weights type
        weights = 'size' ; % 'size' or 'valance'

% Constant objects
    % Boundary nodes
        bndNod = full(mesh.boundaryNodes) ;
        LMBDA = ~bndNod(:).*lmbda(1) + bndNod(:).*lmbda(2) ;
    % Mean over attached elements
        elem2node = mesh.elem2node ;
    % Approximate mesh size
        sz = norm(range(mesh.Nodes,1)) ;
    % Boundaries
        edgeDiff = mesh.edge2node([-1 1]) ;
        bndEdgeDiff = edgeDiff(bndNod,mesh.boundaryEdges) ;
        bndEdge2Node = logical(edgeDiff(bndNod,mesh.boundaryEdges)) ;

% LOOP
    it = 0 ;
    dp2 = Inf ;
    while it<iterations && sqrt(max(dp2))/sz>tol
        dP = zeros([size(mesh.Nodes) 2]) ;
        % Forward/Backward smoothing
            for ii = 1:2 
                % Get the node position
                    X = mesh.Nodes ;
                % Circumcenters
                    Pc = mesh.centroid() ;
                % Boundary edge tangents
                    Te = bndEdgeDiff'*X(bndNod,:) ;
                    Te = Te./sqrt(sum(Te.^2,2)) ;
                    Te = bndEdge2Node*Te ;
                % Normals
                    switch mesh.nCoord
                        case 1
                        case 2
                            Ne = Te*[0 -1 ; 1 0] ;
                        case 3
                            Ne = Te*[0 -1 0 ; 1 0 0 ; 0 0 1] ;
                    end
                    Ne = Ne./sqrt(sum(Ne.^2,2)) ;
                % Weights
                    switch weights
                        case 'valance'
                            W = ones(1,mesh.nElems) ;
                        case 'size' % Areas/Volumes
                            W = mesh.elemSize ;
                    end
                    W = elem2node*sparse(1:mesh.nElems,1:mesh.nElems,W) ;
                % Target point
                    Pt = (W*Pc)./sum(W,2) ; 
                % Updating motion
                    dP(:,:,ii) = Pt-X ;
                    dP(:,:,ii) = dP(:,:,ii).*LMBDA ;
                % Prevent motion of boundary nodes along the normal
                    dP(bndNod,:,ii) = dP(bndNod,:,ii)-sum(dP(bndNod,:,ii).*Ne,2).*Ne ;
                % Update the node position
                    mesh.Nodes = X + sign(1.5-ii) * dP(:,:,ii) ;
            end
        % Norm of displacement
            dp2 = sum(diff(dP,1,3).^2,2) ;
            it = it+1 ;
    end
    

end

