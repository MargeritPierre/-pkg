function mesh = fromBoundary(Bp)
%FROMBOUNDARY Create a simplex mesh from a given list of boundary points

% Compute edge lengths
    dist = sqrt(sum((Bp-circshift(Bp,1,1)).^2,2)) ;
    tooClose = dist<eps*1000*max(dist) ; 
    
% Cull points too close
    Bp(tooClose,:) = [] ;
    dist(tooClose,:) = [] ;

% Parameters
    fixPts = 1:size(Bp,1) ; % fixed points on the boundary
    isInside = @(p)inpolygon(p(:,1),p(:,2),Bp(:,1),Bp(:,2)) ; % tells if a point lies inside the bnd
    targetDist = min(dist) ; % target distance between nodes
    dens = 10./targetDist ; % initial randomized density
    
% Initial distribution of points
    nCoord = size(Bp,2) ;
    bbox = range(Bp,1) ; 
    area = prod(bbox) ;
    N = round(area*dens) ;
    if 0 % random point distribution
        P = rand(N,nCoord) ;
    elseif 1 % regular grid
        N = round((bbox.^nCoord/area*N).^(1/nCoord)) ;
        x = arrayfun(@linspace,N*0,N*0+1,N,'UniformOutput',false) ;
        [x{:}] = ndgrid(x{:}) ;
        x = cat(numel(x)+1,x{:}) ;
        P = reshape(x,[],numel(N)) ;
    end
    P = P.*range(Bp,1) + min(Bp,[],1) ;
    
% Remove outside points
    pts = P(isInside(P),:) ;
    
% Cull point density
    allPts = [Bp ; pts] ;
    % close nodes indices
        [pts,~,ic] = uniquetol(allPts,targetDist,'ByRows',true,'DataScale',1) ;
    % Sparse reprezentation
        M = sparse(ic(:)',1:size(allPts,1),1) ;
    % Force contour points to be fixed
        M(any(M(:,fixPts),2),:) = [] ;
        M = [M ; sparse(1:numel(fixPts),fixPts,1,numel(fixPts),size(M,2))] ;
    % Mean over points
        pts = (M*allPts)./sum(M,2) ;
% Delaunay mesh
    tri = delaunay(pts) ;
    mesh = pkg.geometry.mesh.Mesh('Nodes',pts,'Elems',tri) ;
    mesh.Elems = mesh.Elems.subpart(isInside(mesh.centroid)) ;
    
% Smoothing
    mesh.LaplacianSmoothing([1 0]*.8,10) ;
    
end

