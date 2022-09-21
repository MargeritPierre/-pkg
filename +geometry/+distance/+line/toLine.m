function [d,p,x] = toLine(L1,L2) 
% Compute the shortest distance between infinite lines :
%       between
%   - L1 = cat(3,A1,B1) : [nL nCoord 2]
%       and
%   - L2 = cat(3,A2,B2) : [nL nCoord 2]
% outputs:
%   - line distance d : [nL 1]
%   - parameters p = [p1 p2] : [nL 2]
%   - line points x = cat(3,x1,x2) : [nL nCoord 2]

% We search the minimum distance d between the lines
% We have the line vectors ti = Bi-Ai
% We can parametrize xi(pi) = Ai+ti.pi, pi \in [0 1] , xi(0) = Ai, xi(1) = Bi
% The distance vector D = x2-x1 is such that D.ti = 0 (orthogonal to both lines)
% Which means that D*norm(vectprod(t1,t2)) = d*vectprod(t1,t2)
% (*) D = x2-x1 = (A2-A1) + t2*p2 - t1*p1 = t3/norm(t3)*d with t3 = vectprod(t1,t2)
%   (*).ti : (A2-A1).ti + t2.ti*p2 - t1.ti*p1 = 0 because t3.ti = 0 -> may give pi
%   (*).t3 : (A2-A1).t3  = norm(t3)*d because ti.t3 = 0 -> gives d
%   vectprod((*),t1).t3 : vectprod((A2-A1),t1).t3 - norm(t3)^2*p2 = 0 -> gives p2
%   vectprod((*),t2).t3 : vectprod((A2-A1),t2).t3 - norm(t3)^2*p1 = 0 -> gives p1

if nargin<2 || isempty(L2) ; L2 = L1 ; end

nCoord = size(L1,2) ;
L1 = padarray(L1,[0 3-nCoord],0,'post') ;
L2 = padarray(L2,[0 3-nCoord],0,'post') ;

% Vectors
dA = L2(:,:,1)-L1(:,:,1) ;
t1 = diff(L1,1,3) ; % [nL nCoord]
t2 = diff(L2,1,3) ; % [nL nCoord]
t3 = pkg.math.vectprod(t1,t2) ; % [nL nCoord]
norm3 = sqrt(sum(t3.^2,2)) ; % [nL 1]

% Distance between infinite lines
d = abs(sum(dA.*t3,2))./norm3 ;

% Parameters
if nargout<=1 ; return ; end
p1 = sum(pkg.math.vectprod(dA,t2).*t3,2)./(norm3.^2) ; % [nL 1]
p2 = sum(pkg.math.vectprod(dA,t1).*t3,2)./(norm3.^2) ; % [nL 1]
p = [p1 p2] ; % [nL 2]

% Points
if nargout<=2 ; return ; end
x1 = L1(:,:,1) + t1.*p(:,1) ; % [nL nCoord]
x2 = L2(:,:,1) + t2.*p(:,2) ; % [nL nCoord]
x = cat(3,x1,x2) ; % [nL nCoord 2] ;
x = x(:,1:nCoord,:) ;

end

%% UNITARY TESTS
function tests
%% Random lines
nL = 2 ; nCoord = 2 ;
L1 = rand(nL,nCoord,2) ;
L2 = rand(nL,nCoord,2) ;

[d,p,x] = pkg.geometry.distance.line.toLine(L1,L2) 

L = cat(4,L1,L2) ; % [nL nCoord 2 2]
L = padarray(L,[0 3-nCoord],0,'post') ; % [nL nCoord==3 2 2]
A = permute(L(:,:,1,:),[1 4 2 3]) ; % [nL 2 nCoord==3]
B = permute(L(:,:,2,:),[1 4 2 3]) ; % [nL 2 nCoord==3]
ext = 1/10 ;
Aext = A + (B-A).*(min(p,0)-ext) ; % [nL 2 nCoord==3]
Bext = A + (B-A).*(max(p,1)+ext) ; % [nL 2 nCoord==3]
Lext = reshape(permute(cat(4,Aext,Bext),[4 2 1 3]),2,2*nL,3) ; % [nL*2*2 nCoord==3]
x3 = padarray(x,[0 3-nCoord],0,'post') ;

cla ; axis equal
plot3(A(:,:,1),A(:,:,2),A(:,:,3),'ok')
plot3(B(:,:,1),B(:,:,2),B(:,:,3),'ok')
plot3(Lext(:,:,1),Lext(:,:,2),Lext(:,:,3),'-k')
plot3(permute(x3(:,1,:),[3 1 2]),permute(x3(:,2,:),[3 1 2]),permute(x3(:,3,:),[3 1 2]),'.:r','markersize',25) ;




end