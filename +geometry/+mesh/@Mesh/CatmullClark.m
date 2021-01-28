function mesh = CatmullClark(this,iter)
% Catmull-Clark SURFACE (!) mesh subdivision
% see https://en.wikipedia.org/wiki/Catmull%E2%80%93Clark_subdivision_surface
% ASSUMES THAT mesh.Elems==mesh.Faces !!!!
if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end

if nargin<2 ; iter = 1 ; end

quadType = pkg.geometry.mesh.elements.base.Quadrangle ;

% REPEAT THE SUBDIVISION
for it = 1:iter
    
    % Faces centroids
        Fc = mesh.centroid(mesh.Elems) ;
    % Edges centroids
        Ec = mesh.centroid(mesh.Edges) ;
    % Add the connected face centroids
        Ep = Ec ; % + mesh.elem2edge('mean')*Fc ;
    % New original nodes position
        % Average of each connected face centroid
%             Fn = mesh.elem2node('mean')*Fc ;
        % Average of each connected edge centroid
%             En = (mesh.edge2node*Ec)./N ;
        % New point
            Pn = mesh.Nodes ; %(Fn+2*En+(N-3).*this.Nodes)./N ;
    % Indices of edges in each face
        outEdges = pkg.data.sparse2list(mesh.ElemEdges') ;
        inEdges = outEdges(:,[end,1:end-1]) ;
        nEdgesInElem = sum(outEdges>0,2) ;
        inEdges(:,1) = outEdges(sub2ind(size(outEdges),1:size(outEdges,1),nEdgesInElem(:)')) ;
        cornerNod = mesh.Elems.NodeIdx(:,1:size(outEdges,2)) ;
    % Set the new mesh Nodes
        mesh.Nodes = cat(1,Pn,Ep,Fc) ;
    % Create the quad element indices: 
    %   [cornernode outgoingedgecentroid facecentroid outgoingedgecentroid]
        nodeIdx = cat(3,cornerNod,outEdges,repmat(1:mesh.nElems,[size(cornerNod,2) 1])',inEdges) ;
        nodeIdx = reshape(permute(nodeIdx,[2 1 3]),[],4) ; % Keep approximate element ordering
        nodeIdx(any(isnan(nodeIdx) | nodeIdx==0,2),:) = [] ;
        nodeIdx = double(nodeIdx) + [0 size(Pn,1) size(Pn,1)+size(Ep,1) size(Pn,1)] ;
    % Set the new elements (will clear the connectivity for next iteration)
        mesh.Elems = pkg.geometry.mesh.elements.ElementTable('Types',quadType,'Indices',padarray(nodeIdx,[0 1],1,'pre')) ;
        
end
        

end

