function mesh = cylinder(l_r,dx)
% Surfacic mesh reprezenting a unit cylinder
% l_r can be set as the length to radius ratio

if nargin==0 ; l_r = 1 ; end
if nargin<2 ; dx = .1 ; end

face = pkg.geometry.mesh.shapes.disk(dx) ;
face.nCoord = 3 ;

edge = pkg.geometry.mesh.Mesh('Nodes',face.Nodes,'Elems',face.Edges.subpart(face.boundaryEdges)).clean() ;

if 0 % straight extrusion
    skin = edge.extrude([0 0 l_r],round(l_r/dx)).simplex ;
else % alternate extrusion, to form perfect triangles
    Ne = round(2*l_r/dx/sqrt(3)) ;
% Replicate the edge contour
    x = edge.Nodes + reshape(linspace(0,l_r,Ne+1),1,1,[]).*[0 0 1] ;
% Rotate with half edge length very two replicates
    Le = median(edge.elemSize) ;
    phi = Le/2 ;
    rotMat = pkg.math.rotmat(phi,[0 0 1]) ;
    x(:,:,2:2:end) = permute(pkg.math.innerprod(x(:,:,2:2:end),rotMat,2,2),[1 3 2]) ;
% Generate the triangulation
    edg = double(edge.Elems.NodeIdx) ;
    tri1 = [edg edg(:,1)+edge.nNodes ; flip(edg,2)+edge.nNodes edg(:,2)] + edge.nNodes*reshape(0:2:Ne-1,1,1,[]) ;
    tri2 = [edg edg(:,2)+edge.nNodes ; flip(edg,2)+edge.nNodes edg(:,1)] + edge.nNodes*reshape(1:2:Ne-1,1,1,[])  ;
    tri = cat(3,tri1,tri2) ;
    x = reshape(permute(x,[1 3 2]),[],3) ;
    tri = reshape(permute(tri,[1 3 2]),[],3) ;
    skin = pkg.geometry.mesh.Mesh('Nodes',x,'Elems',tri) ;
end

topface = face.move([0 0 l_r]) ;

face.Elems.NodeIdx = flip(face.Elems.NodeIdx,2) ;

mesh = face + topface + skin ;
mesh = mesh.clean() ;

end