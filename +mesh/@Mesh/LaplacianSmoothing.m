function mesh = LaplacianSmoothing(this,lmbda,iterations,tol)
%Apply Laplacian smoothing to a pkg.mesh.Mesh
% lmbda is a ratio btw 0 and 1; can be [lin lout], to give a different
% value for interior and boundary nodes
if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end

% Format inputs
    % Lambda
        if nargin<2 ; lmbda = [1 0] ; end
        lmbda = lmbda(:)'.*[1 1] ;
    % Iterations
        if nargin<3 ; iterations = Inf ; end
        if nargin<4 ; tol = 1e-9 ; end

% Constant objects
    % Boundary nodes
        bndNod = full(mesh.boundaryNodes) ;
        LMBDA = ~bndNod(:).*lmbda(1) + bndNod(:).*lmbda(2) ;
    % Mean over attached elements
        elem2node = mesh.elem2node ;
        valance = full(sum(elem2node,2)) ;
    % Approximate mesh size
        sz = norm(range(mesh.X.Values,1)) ;

% LOOP
    it = 0 ;
    dp = Inf ;
    while it<iterations && dp/sz>tol
        % Forward Smoothing
            Pc = mesh.circumCenters() ;
            Pt = (elem2node*Pc)./valance ;
            dPf = Pt-mesh.X.Values ;
            dPf = dPf.*LMBDA ;
            mesh.X.Values = mesh.X.Values + dPf ;
        % Backward Smoothing
            Pc = mesh.circumCenters() ;
            Pt = (elem2node*Pc)./valance ;
            dPb = Pt-mesh.X.Values ;
            dPb = dPb.*LMBDA ;
            mesh.X.Values = mesh.X.Values - dPb ;
        % Norm of displacement
            dp = max(sqrt(sum((dPf-dPb).^2,2))) ;
            it = it+1 ;
    end
    

end

