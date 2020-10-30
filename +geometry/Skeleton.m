classdef Skeleton < pkg.geometry.mesh.Mesh
%SKELETON the skeleton of a geometry

properties %(SetAccess = protected)
    BoundaryPoints % [nPts nCoord] list of points defining the object boundaries
    InsideFcn % function handle saying if a points lies inside the contour
    DelaunayMesh % Delaunay triangulation generated for the skeleton
end

methods
    
    function this = Skeleton(input)
    % Object constructor
        if nargin==0 ; return ; end
        if isa(input,'pkg.geometry.mesh.Mesh')
            this.BoundaryPoints = input.Nodes(input.BoundaryNodes,:) ;
            this.InsideFcn = @input.isInside ;
        elseif isa(input,'pkg.geometry.levelset.LevelSet')
            this.BoundaryPoints = input.discretizeContour ;
            this.InsideFcn = @input.inside ;
        else % Arbitrary points
            this.BoundaryPoints = input ;
            this.InsideFcn = @(P)inpolygon(P(:,1),P(:,2),this.BoundaryPoints(:,1),this.BoundaryPoints(:,2)) ;
        end
        this.buildSkeleton ;
    end
    
    function buildSkeleton(this)
    % Build the skeleton from BoundaryPoints and Inside function
    % Create the delaunay triangulation
        tri = delaunay(this.BoundaryPoints) ;
        elmtType = pkg.geometry.mesh.elements.base.Triangle ;
        indices = padarray(tri,[0 1],1,'pre') ;
        elems = pkg.geometry.mesh.elements.ElementTable('Types',elmtType,'Indices',indices) ;
        this.DelaunayMesh = pkg.geometry.mesh.Mesh('Nodes',this.BoundaryPoints,'Elems',elems) ;
    % Delaunay mesh cleaning
        % Remove triangles that lie outside the object
            valid = this.InsideFcn(this.DelaunayMesh.centroid) ;
        % Remove triangles with a too small area
            X = this.DelaunayMesh.Elems.dataAtIndices(this.DelaunayMesh.Nodes) ;
            areas = polyarea(X(:,:,1),X(:,:,2),2) ;
            valid = valid & areas > 1e-6 * median(areas(:)) ;
        % Keep only valid triangles
            this.DelaunayMesh.Elems = this.DelaunayMesh.Elems.subpart(valid) ;
    % Compute element circumcenters
        X = this.DelaunayMesh.Elems.dataAtIndices(this.DelaunayMesh.Nodes) ;
        this.Nodes = pkg.math.circumcenter(permute(X,[2 3 1])) ;
    % Create edges linking circumcenters 
    % (one link by inner edge in the delaunay mesh)
        elem2edge = this.DelaunayMesh.elem2edge ;
        innerEdg = sum(elem2edge,2)==2 ;
        links = pkg.data.sparse2list(elem2edge(innerEdg,:)) ;
        elmtType = pkg.geometry.mesh.elements.base.Bar ;
        indices = padarray(links,[0 1],1,'pre') ;
        this.Elems = pkg.geometry.mesh.elements.ElementTable('Types',elmtType,'Indices',indices) ;
    end
end


methods
% In all the following, indices are as follows:
%   - idx <= this.nNodes denote a node of the skeleton (==delaunay elem circum.)
%   - idx > this.nNodes denote a node on the contour (==delauney mesh node)
    
    function X = allNodes(this)
    % Return all the node coordinates corresponding to given indices
        X = [this.Nodes ; this.DelaunayMesh.Nodes] ;
    end
    
    function idx = bndEdges(this)
    % Return boundary edges
        outEdg = this.DelaunayMesh.Edges.subpart(this.DelaunayMesh.outerEdges) ;
        idx = outEdg.NodeIdx + this.nNodes ;
    end    
    
    function idx = skelEdges(this)
    % Create a graph representation of the skeleton
        idx = this.Edges.NodeIdx ;
    end 
    
    function idx = singularIdx(this)
    % Indices of singular elements or corresponding circumcenters
    % Singular triangles are linked to only inner edges
        elem2edge = this.DelaunayMesh.elem2edge ;
        innerEdg = sum(elem2edge,2)==2 ;
        idx = find((innerEdg'*elem2edge)==3) ;
    end
    
    function idx = cornerIdx(this)
    % Indices of corner nodes on the boundary
    % Corner nodes are linked to boundary edges only
        idx = find(this.DelaunayMesh.endNodes) + this.nNodes ;
    end
        
    function mesh = quadPatches(this)
    % Return a collection of quad patch representing the skeleton object
    % skeleton & boundary topology
    % see R. Oval & al., Topology Finding of Structural Patterns
    % https://hal.archives-ouvertes.fr/hal-01883505/document
    % node indices > this.nNodes represent countourpoint indices
        X = [this.Nodes ; this.DelaunayMesh.Nodes] ;
    % S-points, circumcenters of singular triangles
    % (that are related to only inner edges)
        elem2edge = this.DelaunayMesh.elem2edge ;
        innerEdg = sum(elem2edge,2)==2 ;
        S = (innerEdg'*elem2edge)==3 ;
        Sidx = find(S) ;
    % B-points, as the vertices of the singular triangles
        Bidx = this.DelaunayMesh.Elems.NodeIdx(S,:) + this.nNodes ;
        linkSB = reshape(cat(3,repmat(Sidx(:),[1 3]),Bidx),[],2) ;
    % Singular contour points, linked to only one element
        C = this.DelaunayMesh.endNodes ;
    % Associated triangles
        [Cidx,Ctri] = find(C(:) & this.DelaunayMesh.elem2node) ;
        linkC = [Cidx(:)+this.nNodes Ctri(:)] ;
    % Return the mesh
        newIdx = [linkSB ; linkC] ;
        indices = padarray([this.Elems.NodeIdx(:,1:2) ; newIdx],[0 1],1,'pre') ;
        elmtType = pkg.geometry.mesh.elements.base.Bar ;
        elems = pkg.geometry.mesh.elements.ElementTable('Types',elmtType,'Indices',indices) ;
        mesh = pkg.geometry.mesh.Mesh('Nodes',X,'Elems',elems) ;
    end
end


end

