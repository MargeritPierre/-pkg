function X = mtimes(A,B)
%MTIMES Return the products X of stacks of matrices A and B
% Applies element-wise operations if possible
% 
%   X = MTIMES(A,B)
%       where size(A) = [n k sz1 sz2 sz3 ..]
%             size(B) = [k m sz1 sz2 sz3 ..] !


% Process inputs
    szA = size(A) ; 
    szB = size(B) ; 
    [n,k,m] = deal(szA(1),szA(2),szB(2)) ;
    if szB(1)~=k ; error('Incompatible sizes for A:[n k] and B:[k m]') ; end
    nDim = max(ndims(A),ndims(B)) ;
    SZ = ones(1,nDim+2) ;
    SZ(1:ndims(A)) = szA ;
    SZ(1:ndims(B)) = max(SZ(1:ndims(B)),szB) ;
    SZ(1:2) = [] ;
    A = permute(A,[1 nDim+1 3:nDim 2]) ; % [n 1 sz1 sz2 ... k]
    B = permute(B,[nDim+1 2 3:nDim 1]) ; % [1 m sz1 sz2 ... k]

% Compute the product by permutation-summation
    X = sum(A.*B,ndims(A)) ;
    
% Reshape X
    X = reshape(X,[n m SZ 1 1]) ;

return


%% TESTS

%% ONE INPUT iA = mtimes(A)
    n = 3 ; k = 3 ; m = 3 ; N = 10000 ;
    A = rand(n,k,N)  + 1i*rand(n,k,N) ;
    B = rand(k,m,N)  + 1i*rand(k,m,N) ;
    tic ; X0 = pkg.data.matrixfun(@mtimes,A,B) ; t0 = toc ;
    tic ; X = pkg.math.mtimes(A,B) ; t = toc 
    disp(['PINV(A): ' ...
            's=[' num2str(n) ',' num2str(m) ']' ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(X(:)-X0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ])

