%% COMPUTATIONNAL EFFICIENCY OF SYSTEM MATRIX ASSEMBLY
clc ; clearvars
P = rand(100000,3) ;
tri = delaunay(P) ;

mesh = pkg.geometry.mesh.Mesh('Nodes',P,'Elems',tri) ;

[ee,we,ie] = mesh.integration();
N = mesh.interpMat(ee,ie) ;
W = diag(sparse(we)) ;

A = N' ; B = W ; C = N ;

% Brute force
tic ; K = A*B*C ; toc

% % Kronecker: vec(A*B*C) = kron(C',A)*vec(B)
% tic ; 
%     CxA1 = kron(C',A) ; 
%     k1 = CxA1*B(:) ; 
% toc
% norm(K(:)-k1)

% If W is diagonal, vec(A*diag(b)*C) = khatrirao(C',A)*b
tic ; 
    b = diag(B) ;
    CxA2 = pkg.math.khatrirao(C',A) ;
    k2 = CxA2*b ; 
toc
% norm(K(:)-k2)
