function X = mldivide(A,B)
%MLDIVIDE Return the solutions X of stacks of small matrices A and B
% Applies element-wise operations if possible
% 
%   X = MLDIVIDE(A,B)
%       where size(A) = [n p sz1 sz2 sz3 ..]
%             size(B) = [n k sz1 sz2 sz3 ..] !

% Case where A matrix and B vector (k==1)
    if ndims(B)==ndims(A)-1
        szB = size(B) ;
        B = reshape(B,[szB(1) 1 szB(2:end)]) ;
    end

% Process inputs
    szA = size(A) ; 
    szB = size(B) ;
    [n,p,k] = deal(szA(1),szA(2),szB(2)) ;
    if szB(1)~=n ; error('Incompatible sizes for A:[n p] and B:[n k]') ; end
    SZ = [szA 1 1] ; SZ(1:2) = [] ;
    A = A(:,:,:) ; B = B(:,:,:) ;

% Solve A.X = B with the pseudo-inverse of A
    piA = pkg.math.pinv(A) ;
    X = pkg.math.mtimes(piA,B) ;
    
% Reshape X
    X = reshape(X,[p k SZ]) ;

return


%% TESTS

%% ONE INPUT iA = mldivide(A)
    n = 3 ; p = 3 ; k = 1 ; N = 10000 ;
    A = rand(n,p,N)  + 1i*rand(n,p,N) ;
    B = rand(n,k,N)  + 1i*rand(n,k,N) ;
    tic ; X0 = pkg.data.matrixfun(@mldivide,A,B) ; t0 = toc ;
    tic ; X = pkg.math.mldivide(A,B) ; t = toc 
    disp(['PINV(A): ' ...
            's=[' num2str(n) ',' num2str(m) ']' ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(X(:)-X0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ])

