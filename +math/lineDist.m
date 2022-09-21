function [d,p,x] = lineDist(L1,L2,isinf) 
% Compute the shortest distance between infinite lines (or finite segments) :
%       between
%   - L1 = cat(3,A1,B1) : [nL nCoord 2]
%       and
%   - L2 = cat(3,A2,B2) : [nL nCoord 2]
% isinf: boolean, false for finite segments
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
if nargin<3 ; isinf = false ; end

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
if isinf ; d = abs(sum(dA.*t3,2))./norm3 ; end

% Parameters
if isinf && nargout<=1 ; return ; end
p1 = sum(pkg.math.vectprod(dA,t2).*t3,2)./(norm3.^2) ; % [nL 1]
p2 = sum(pkg.math.vectprod(dA,t1).*t3,2)./(norm3.^2) ; % [nL 1]
p = [p1 p2] ; % [nL 2]

% Finite segments
if ~isinf 
    isout = any(p<0 | p>1,2) ; 
    
end

% Points
if isinf && nargout<=2 ; return ; end
x1 = L1(:,:,1) + t1.*p(:,1) ; % [nL nCoord]
x2 = L2(:,:,1) + t2.*p(:,2) ; % [nL nCoord]
x = cat(3,x1,x2) ; % [nL nCoord 2] ;
x = x(:,1:nCoord,:) ;

if ~isinf ; d = sqrt(sum(diff(x,1,3).^2,2)) ; end

end

%% UNITARY TESTS
function tests
%% Random lines
nL = 10 ; nCoord = 3 ;
L1 = rand(nL,nCoord,2) ;
L2 = rand(nL,nCoord,2) ;

[d,p,x] = pkg.math.lineDist(L1,L2) 

cla ; axis equal
plot3(permute(L1(:,1,:),[3 1 2]),permute(L1(:,2,:),[3 1 2]),permute(L1(:,3,:),[3 1 2]),'.-b','markersize',25) ;
plot3(permute(L2(:,1,:),[3 1 2]),permute(L2(:,2,:),[3 1 2]),permute(L2(:,3,:),[3 1 2]),'.-r','markersize',25) ;
plot3(permute(x(:,1,:),[3 1 2]),permute(x(:,2,:),[3 1 2]),permute(x(:,3,:),[3 1 2]),'.:k','markersize',25) ;




end