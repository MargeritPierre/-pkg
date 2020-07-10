function mesh = CatmullClark(this,iter)
% Catmull-Clark SURFACE mesh subdivision
% see https://en.wikipedia.org/wiki/Catmull%E2%80%93Clark_subdivision_surface
if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end

if nargin<2 ; iter = 1 ; end

% REPEAT THE SUBDIVISION
for it = 1:iter
    
    % Faces centroids
        Fc = mesh.centroid(mesh.Faces) ;
    % Edges centroids
        Ec = mesh.centroid(mesh.Edges) ;
    % Add the connected face centroids
        face2edge = mesh.face2edge('mean') ;
        Ep = Ec ; %+ face2edge*Fc ;
    % New original nodes position
        % Average of each connected face centroid
            face2node = mesh.face2node('mean') ;
            N = sum(face2node,2) ;
            Fn = (face2node*Fc)./N ;
        % Average of each connected edge centroid
            En = (mesh.edge2node*Ec)./N ;
        % New point
            Pn = mesh.X.Values ; %(Fn+2*En+(N-3).*this.X.Values)./N ;
    % Indices of edges in each face
        [~,FaceEdges] = mesh.Edges.NodeIdx ;
    % Build the new mesh
        elems = NaN(0,4) ;
        for pp = 1:size(mesh.Elems,2)
            if pp==size(mesh.Elems,2) 
                %lastInd = FaceEdges(sub2ind(size(FaceEdges),(1:this.nElems)',sum(~isnan(this.Elems),2))) ;
                lastInd = FaceEdges(sub2ind(size(FaceEdges),(1:mesh.nElems)',ones(mesh.nElems,1))) ; 
            else 
                lastInd = FaceEdges(:,pp+1) ; 
            end
            elems = [ elems ; ...
                        mesh.Elems(:,pp) ...
                        FaceEdges(:,pp)+mesh.nNodes ...
                        (1:mesh.nElems)'+mesh.nNodes+mesh.nEdges ...
                        lastInd(:)+mesh.nNodes ...
                        ] ;
        end
    % Cull non-quads elems
        elems(any(isnan(elems),2),:) = [] ;
    % Set the mesh
        mesh.Nodes = cat(1,Pn,Ep,Fc) ;
        mesh.Elems = elems ;
        
end
        

end

