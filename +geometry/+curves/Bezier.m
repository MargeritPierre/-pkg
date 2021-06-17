%% Compute Bezier Curve
% see https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-construct.html

% Parameters
n = 3 ; % curve degree
m = n+1 ; % number of control points
u = linspace(0,1,1000)' ; % evaluation abscissa

% Cumpute the function basis
ii = 0:n ;
Bni = (factorial(n)./(factorial(ii).*factorial(n-ii))).*(u.^ii).*((1-u).^(n-ii)) ;
clf ; plot(u,Bni) ;

% Evaluate a curve
C = rand(m,2) ;
x = Bni*C ;
clf ; 
plot(x(:,1),x(:,2),'r') ; 
plot(C(:,1),C(:,2),'.-.k','linewidth',1.5,'markersize',20)
