%% TESTS ON ITERATIVE SOLVERS

%% STEADY TEMPERATURE TEST
L = [1 1] ;
x0 = [.5 .5].*L ;
Q = 1 ;
lambda = 1 ;
dx0 = norm(L)/20 ;
subdiv = 1 ;

lvlst = pkg.geometry.levelset.Rectangle(L.*[0;1]) ;
lvlst.Kinks(end+1,:) = x0 ;
mesh = lvlst.mesh(dx0) ;
mesh.sortElems ;
for ss = 1:subdiv
    nodIdx = [mesh.Elems.NodeIdx mesh.nNodes + pkg.data.sparse2list(mesh.ElemEdges')] ;
    nodIdx = reshape(nodIdx(:,[1 4 6 ; 2 5 4 ; 3 6 5 ; 4 5 6]),[],3) ;
    %nodIdx = reshape(nodIdx(:,[4 5 6]),[],3) ;
    mesh.Nodes = [mesh.Nodes ; mesh.centroid('Edges')] ;
    mesh.Elems.Indices = padarray(nodIdx,[0 1],1,'pre') ;
end
%mesh.CatmullClark(1) ;

clf ; axis equal tight
pl = plot(mesh) ;
%plot(x0(1),x0(2),'or')

[ee,w,ie] = mesh.integration() ;
W = diag(sparse(w)) ;
G = mesh.gradMat(1,ee,ie) ;

keepDOF = ~mesh.near([0 NaN]) ;

K = G'*kron(speye(2)*lambda,W)*G ;
%K = .5*(K+K') ;
f = sparse(pkg.data.closestPoint(x0,mesh.Nodes),1,Q,mesh.nNodes,1) ;

T = zeros(mesh.nNodes,1) ;
A = K(keepDOF,keepDOF) ;
b = f(keepDOF) ;


%%
tol = 1e-10 ;
maxIt = 1000 ;
options = struct('type','nofill','michol','on','diagcomp',.1) ;

tic ; 
    t0 = A\b ; 
toc

tic ; 
%     [PP,RR,CC] = equilibrate(A) ;
%     A = RR*PP*A*CC ;
%     b = RR*PP*b ;
    L = ichol(A,options) ;
    t1 = pcg(A,b,tol,maxIt,L,L') ; 
toc
err1 = norm(t0-t1)/norm(t0)


T(keepDOF) = t0 ;

pl.CData = T ;


