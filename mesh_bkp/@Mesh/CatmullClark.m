function this = CatmullClark(this,iter)
% Catmull-Clark SURFACE mesh subdivision
% see https://en.wikipedia.org/wiki/Catmull%E2%80%93Clark_subdivision_surface

if nargin<2 ; iter = 1 ; end

% REPEAT THE SUBDIVISION
for it = 1:iter
    
    % Faces circumcenters
        Fc = this.circumCenters(this.Faces) ;
    % Edges circumcenters
        Ec = this.circumCenters(this.Edges) ;
    % Add the connected face circumcenters
        face2edge = this.face2edge('mean') ;
        Ep = Ec ; %+ face2edge*Fc ;
    % New original nodes position
        % Average of each connected face circumcenter
            face2nod = this.face2nod('mean') ;
            N = sum(face2nod,2) ;
            Fn = (face2nod*Fc)./N ;
        % Average of each connected edge circumcenter
            En = (this.edge2nod*Ec)./N ;
        % New point
            Pn = this.Nodes ; %(Fn+2*En+(N-3).*this.Nodes)./N ;
    % Indices of edges in each face
        [~,FaceEdges] = this.listEdges(this.Elems) ;
    % Build the new mesh
        elems = NaN(0,4) ;
        for pp = 1:size(this.Elems,2)
            if pp==size(this.Elems,2) 
                %lastInd = FaceEdges(sub2ind(size(FaceEdges),(1:this.nElems)',sum(~isnan(this.Elems),2))) ;
                lastInd = FaceEdges(sub2ind(size(FaceEdges),(1:this.nElems)',ones(this.nElems,1))) ; 
            else 
                lastInd = FaceEdges(:,pp+1) ; 
            end
            elems = [ elems ; ...
                        this.Elems(:,pp) ...
                        FaceEdges(:,pp)+this.nNodes ...
                        (1:this.nElems)'+this.nNodes+this.nEdges ...
                        lastInd(:)+this.nNodes ...
                        ] ;
        end
    % Cull non-quads elems
        elems(any(isnan(elems),2),:) = [] ;
    % Set the mesh
        this.Nodes = cat(1,Pn,Ep,Fc) ;
        this.Elems = elems ;
        
end
        

end

