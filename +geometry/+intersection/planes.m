function [p,ind] = planes(O,N)
%PLANES Intersection points between multiple planes
%   - O [nPl nCoord] % plane origin(s)
%   - N [nPl nCoord] % plane normal(s)
% outputs: 
%   - p [nComb nCoord] % intersection points
%   - ind [nComb nCoord] % indices of corresponding planes
% nComb is the number of intersection combinations (see below)

[nPl,nCoord] = size(O) ;
if nPl<nCoord ; error(['At least ' num2str(nCoord) ' planes must be given to determine a ' num2str(nCoord) 'D intersection point']) ; end

% The intersection of nCoord planes defines a point in (nCoord)D-space
% Get unique combination of (non-repeated) planes
ind = pkg.math.combinations(repmat({1:nPl},[1 nCoord])) ;
ind = unique(sort(ind,2),'rows') ; % unique combinations
ind(any(diff(ind,1,2)==0,2),:) = [] ; % non-repeating planes
nComb = size(ind,1) ;

% Determine intersections
% A point p is on a plane if dot(p-O,N)=0 => N*p' = N*O' ; 
O = permute(reshape(O(ind,:),[nComb nCoord nCoord]),[2 3 1]) ;
N = permute(reshape(N(ind,:),[nComb nCoord nCoord]),[2 3 1]) ;
p = NaN(nComb,nCoord) ;
for cc = 1:nComb
    n = N(:,:,cc) ; o = O(:,:,cc) ;
    p(cc,:) = n\sum(n.*o,2) ;
end

end

function test
%% TEST FUNCTION
nCoord = 2 ; nPl = 4 ;
O = 2*rand(nPl,nCoord)-1 ; [1 0 ; 0 1] ;
N = 2*rand(nPl,nCoord)-1 ; [1 0 ; 1 0] ; 

tic
[p,ind] = pkg.geometry.intersection.planes(O,N) ;
toc

clf ;
axis equal
p3D = padarray(p,[0 max(0,3-nCoord)],0,'post') ;
O3D = padarray(O,[0 max(0,3-nCoord)],0,'post') ;
N3D = padarray(N,[0 max(0,3-nCoord)],0,'post') ;
quiver3(O3D(:,1),O3D(:,2),O3D(:,3),N3D(:,1),N3D(:,2),N3D(:,3),'k') ;
T3D = permute(O3D + 4*cat(3,-1,1).*N3D(:,[2 1 3]).*[-1 1 1],[3 1 2]) ;
plot3(T3D(:,:,1),T3D(:,:,2),T3D(:,:,3),'-.k','linewidth',1)
plot3(p3D(:,1),p3D(:,2),p3D(:,3),'.k','markersize',20)

end

