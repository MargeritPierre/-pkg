%% ADVANCING FRONT MESH
    tag = 'advancingfront' ;
    delete(findobj(gca,'tag',tag)) ;
    mesh = pkg.geometry.mesh.Mesh() ;
    pl = plot(mesh,'Tag',tag) ;
    pl.VisibleNodes = 'all' ;
    pl.Nodes.MarkerSize = 10 ;
    
    fd = @lvlst.Function ;
    fgrad = @lvlst.gradient ;
    
    h0 = norm(range(lvlst.BoundingBox,1))/50 ;
    fh = @(p)repmat(h0,[size(p,1) 1]) ;
    
    % Start from the discretized contour
        mesh.Nodes = lvlst.discretizeContour(fh);
        lastNodes = 1:mesh.nNodes ;
    
    % Edges on the front: initialize
        mesh.Elems.Types = pkg.geometry.mesh.elements.base.Triangle ;
%         mesh.Elems.NodeIdx = delaunay(mesh.Nodes) ;
%         mesh.Elems = mesh.Elems.subpart(lvlst.Function(mesh.centroid)<0) ;
    
    maxIt = 2 ; 2*norm(range(lvlst.BoundingBox,1))/h0 ;
    deps = 1e-3 ;
    it = 0 ;
    while ~isempty(lastNodes) && it<maxIt
        it = it+1 ;
        p = mesh.Nodes(lastNodes,:) ;
    % Keep a copy of the current points
        p0 = p ;
        d0 = fd(p) ;
    % New layer of points
        for iii = 1:10
        % Direction of the gradient
            grad = fgrad(p) ;
        % Target edge length
            ht = fh((p0+p)/2) ;
        % Current edge length
            h = d0-fd(p) ; %sqrt(sum((p0-p).^2,2)) ;
        % Move accordingly
            p = p-grad.*(ht-h) ;
        end
    % Cull points for which the signed distance increased
        d = fd(p) ;
        p(d>d0,:) = [] ;
        d0 = d(d>d0) ;
    % Add Nodes to the mesh
        lastNodes = mesh.nNodes + (1:size(p,1)) ;
        mesh.Nodes = [mesh.Nodes ; p] ;
    % Display
        if 1 ; pl.update ; drawnow ; end
    end
    
    pl.update
    
    