function [d,p,x] = toSegment(S1,S2) 
% Compute the shortest distance between finite segments :
%       between
%   - S1 = cat(3,A1,B1) : [nS nCoord 2]
%       and
%   - S2 = cat(3,A2,B2) : [nS nCoord 2]
% outputs:
%   - line distance d : [nS 1]
%   - parameters p = [p1 p2] : [nS 2]
%   - line points x = cat(3,x1,x2) : [nS nCoord 2]

if nargin<2 || isempty(S2) ; S2 = S1 ; end

nCoord = size(S1,2) ;
S1 = padarray(S1,[0 3-nCoord],0,'post') ;
S2 = padarray(S2,[0 3-nCoord],0,'post') ;

% Inifinite line distance
[d,p,x] = pkg.geometry.distance.line.toLine(S1,S2) ;

% Is the parameter in the authorized domain ?
isout = any( p<0 | p>1 ,2) ;
if ~any(isout) ; return; end
d(isout) = Inf ;

% Distance from one segment ends to the other segment
[da1s2,pa1s2] = pkg.geometry.distance.point.toSegment(S1(:,:,1),S2(:,:,1),S2(:,:,2)) ; % [nL 1]
[db1s2,pb1s2] = pkg.geometry.distance.point.toSegment(S1(:,:,2),S2(:,:,1),S2(:,:,2)) ; % [nL 1]
[da2s1,pa2s1] = pkg.geometry.distance.point.toSegment(S2(:,:,1),S1(:,:,1),S1(:,:,2)) ; % [nL 1]
[db2s1,pb2s1] = pkg.geometry.distance.point.toSegment(S2(:,:,2),S1(:,:,1),S1(:,:,2)) ; % [nL 1]

% Minimum distance
[d,imin] = min([d da1s2 db1s2 da2s1 db2s1],[],2) ;

% Associated parameter
p(imin==2,:) = pa1s2(imin==2).*[0 1] + [0 0] ; 
p(imin==3,:) = pb1s2(imin==3).*[0 1] + [1 0] ; 
p(imin==4,:) = pa2s1(imin==4).*[1 0] + [0 0] ; 
p(imin==5,:) = pb2s1(imin==5).*[1 0] + [0 1] ; 

% Points
x1 = S1(:,:,1) + diff(S1,1,3).*p(:,1) ; % [nL nCoord]
x2 = S2(:,:,1) + diff(S2,1,3).*p(:,2) ; % [nL nCoord]
x = cat(3,x1,x2) ; % [nL nCoord 2] ;
x = x(:,1:nCoord,:) ;

end

%% UNITARY TESTS
function tests
%% Random lines
nL = 200000 ; nCoord = 3 ;
L1 = rand(nL,nCoord,2) ;
L2 = rand(nL,nCoord,2) ;

tic ; [d,p,x] = pkg.geometry.distance.segment.toSegment(L1,L2) ; toc
err = norm(d-sqrt(sum(diff(x,1,3).^2,2)))/norm(d)

%%

L = cat(4,L1,L2) ; % [nL nCoord 2 2]
L = padarray(L,[0 3-nCoord],0,'post') ; % [nL nCoord==3 2 2]
L = reshape(permute(L,[3 4 1 2]),2,2*nL,3) ; % [2 2*nL nCoord==3]
x3 = padarray(x,[0 3-nCoord],0,'post') ;

cla ; axis equal
plot3(L(:,:,1),L(:,:,2),L(:,:,3),'.-k','markersize',25)
plot3(permute(x3(:,1,:),[3 1 2]),permute(x3(:,2,:),[3 1 2]),permute(x3(:,3,:),[3 1 2]),'.:r','markersize',25) ;




end