function mesh = rectangle(sz,dx)
% Surfacic mesh reprezenting a rectangle
if nargin<1 ; sz = [sqrt(2) 1] ; end % default rectangle
if nargin<2 ; dx = sz ; end % by default return very simple meshes

mesh = pkg.geometry.mesh.GridMesh(sz.*[0;1],dx).simplex ;

end

