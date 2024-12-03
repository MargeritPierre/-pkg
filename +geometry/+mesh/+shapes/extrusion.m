function mesh = extrusion(surfmesh,vec,dx)
% Surfacic mesh reprezenting a unit cylinder
% l_r can be set as the length to radius ratio

if nargin==0 ; dx = .1 ; end
if nargin<2 ; l_r = 1 ; end

bottomface = pkg.geometry.mesh.shapes.disk(dx) ;



end