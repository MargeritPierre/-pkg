function this = LaplacianSmoothing(this,lmbda,iterations,tol)
%Apply Laplacian smoothing to a pkg.mesh.Mesh
% lmbda is a ratio btw 0 and 1; can be [lin lout], to give a different
% value for interior and boundary nodes

% Format inputs
    % Lambda
        if nargin<2 ; lmbda = [1 0] ; end
        lmbda = lmbda(:)'.*[1 1] ;
    % Iterations
        if nargin<3 ; iterations = Inf ; end
        if nargin<4 ; tol = 1e-9 ; end

% Constant objects
    % Boundary nodes
        bndNod = full(this.boundaryNodes) ;
        LMBDA = ~bndNod(:).*lmbda(1) + bndNod(:).*lmbda(2) ;
    % Mean over attached elements
        elem2nod = this.elem2nod ;
        valance = full(sum(elem2nod,2)) ;
    % Approximate mesh size
        sz = norm(range(this.Nodes,1)) ;

% LOOP
    it = 0 ;
    dp = Inf ;
    while it<iterations && dp/sz>tol
        % Forward Smoothing
            Pc = this.circumCenters() ;
            Pt = (elem2nod*Pc)./valance ;
            dPf = Pt-this.Nodes ;
            dPf = dPf.*LMBDA ;
            this.Nodes = this.Nodes + dPf ;
        % Backward Smoothing
            Pc = this.circumCenters() ;
            Pt = (elem2nod*Pc)./valance ;
            dPb = Pt-this.Nodes ;
            dPb = dPb.*LMBDA ;
            this.Nodes = this.Nodes - dPb ;
        % Norm of displacement
            dp = max(sqrt(sum((dPf-dPb).^2,2))) ;
            it = it+1 ;
    end
    

end

