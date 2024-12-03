function mesh = box(sz,dx)
% Surfacic mesh reprezenting a rectangle
if nargin<1 ; sz = [sqrt(2) 1 sqrt(3)] ; end % default rectangle
if nargin<2 ; dx = sz ; end % by default return very simple meshes
dx = dx.*[1 1 1] ;

face = pkg.geometry.mesh.shapes.rectangle(sz(1:2),dx(1:2)) ;
face.nCoord = 3 ;

edge = pkg.geometry.mesh.Mesh('Nodes',face.Nodes,'Elems',face.Edges.subpart(face.boundaryEdges)).clean() ;
skin = edge.extrude([0 0 sz(3)],round(sz(3)/dx(3))) ;

topface = face.move([0 0 sz(3)]) ;

face.Elems.NodeIdx = flip(face.Elems.NodeIdx,2) ;

mesh = face + topface + skin ;
mesh = mesh.clean().simplex ;

end

