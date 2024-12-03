function mesh = disk(dx)
% Surfacic mesh reprezenting a unit disk

if nargin==0 ; dx = .1 ; end

Nr = round(2/dx/sqrt(3)) ;

circumrad = dx/sqrt(3) ;
R = linspace(circumrad,1,Nr) ;
Nt = floor(2*pi*R/dx) ;

r = repelem(R,Nt) ;
dt = repelem(2*pi./Nt,Nt) ;
t = cumsum(dt)-repelem(2*pi*(0:Nr-1),Nt) ;

x = r(:).*[cos(t(:)) sin(t(:))] ;
tri = delaunay(x) ;

mesh = pkg.geometry.mesh.Mesh('Nodes',x,'Elems',tri) ;

end