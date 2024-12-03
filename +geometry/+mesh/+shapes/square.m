function mesh = square(dx)
% Surfacic mesh reprezenting a square

if nargin<1 ; dx = 1 ; end % by default return very simple meshes
mesh = pkg.geometry.mesh.shapes.rectangle([1 1],dx) ;

end