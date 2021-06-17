%% Compute BSpline functions with their derivatives
% see https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-coef.html
clc

q = 6 ; % number of control points
p = 3 ; % polynomial degree
der = 0 ; % derivation order

m = q+p+1 ; % number of knots
ui = padarray((0:q-p)/(q-p),[0 p],'replicate','both') ; % knot vector

u = linspace(0,1,100)' ; % evaluation coordinates
u = linspace(-1,2,1001)' ; % evaluation coordinates

% Zero-th order basis (door) function
N = u>=ui(1:end-1) & u<ui(2:end) ;
N(u==1,end-p) = 1 ;
clf ; plot(u,N) ;

% Recursion formula
dN = N ;
for pp = 1:p
    kk = 1:size(N,2)-1 ; % knot spans
    N = (u-ui(kk))./max(ui(kk+pp)-ui(kk),eps).*N(:,1:end-1) ...
        + (ui(kk+pp+1)-u)./max(ui(kk+pp+1)-ui(kk+1),eps).*N(:,2:end) ;
    if pp==p-der ; dN = N ; end
end
clf ; plot(u,N) ;

if der==0 ; return ; end

% Analytical derivative
for pp = p-der+1:p
    kk = 1:size(dN,2)-1 ; % knot spans
    dN = pp./max(ui(kk+pp)-ui(kk),eps).*dN(:,1:end-1) ...
        - pp./max(ui(kk+pp+1)-ui(kk+1),eps).*dN(:,2:end) ;
end
clf ; plot(u,dN) ;

% Finite differences
ddN = N ;
for dd = 1:der
    ddN = [...
            diff(ddN([1 2],:),1,1) ; ...
            (ddN(3:end,:)-ddN(1:end-2,:))/2 ; ...
            diff(ddN(end-1:end,:),1,1) ...
          ]/mean(diff(u)) ;
end
set(gca,'ColorOrderIndex',1) ; plot(u,dN,'.','markersize',15) ;


