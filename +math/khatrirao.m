function M = khatrirao(varargin)
%KHATRIRAO Computes the khatri-rao product of a suite of matrices
% varargin contains Q>=2 matrices U{1}...U{Q}
% U{q} is of size [Nq K] , K is common to every matrix 
% M is of size [prod(Nq) K]
    if numel(varargin)==1 ; M = varargin{1} ; return ; end

% Check sizes
    [~,K] = cellfun(@size,varargin) ;
    if range(K)~=0 ; error('The matrices must share the same number of columns') ; end
    
% Recursive call
    A = varargin{1} ;
    if nargin==2 ; B = varargin{2} ;
    else ; B = pkg.math.khatrirao(varargin{2:end}) ;
    end

% Depending on the sparsity..
    if issparse(A) && issparse(B)
        M = sparseProd(A,B) ;
    else
        M = fullProd(A,B) ;
    end

end

% KHATRI-RAO PRODUCT OF TWO FULL MATRICES
function M = fullProd(A,B)
% KR product of full matrices
    M = reshape(permute(B,[1 3 2]).*permute(A,[3 1 2]),[],size(A,2)) ;
end

% KHATRI-RAO PRODUCT OF TWO SPARSE MATRICES
function M = sparseProd(A,B)
% KR product of full matrices
    [ia,ja,va] = find(A) ;
    [ib,jb,vb] = find(B) ;
    [aa,bb] = pkg.data.findIn(ja,jb) ; 
    jj = ja(aa) ; % all ja(aa)==jb(bb)
    ii = (ia(aa)-1)*size(B,1) + ib(bb) ;
    vv = va(aa).*vb(bb) ;
    M = sparse(ii,jj,vv,size(A,1)*size(B,1),size(A,2)) ;
end



%% UNIT TESTS
function tests
%% FULL MATRICES
clc ; clearvars ;
Nq = [1000 200] ;
K = 100 ;
U = {} ; for mm=1:numel(Nq) ; U{mm} = rand([Nq(mm) K]) ; end
tic ; M1 = pkg.math.khatrirao(U{:}) ; toc
tic ; M2 = repelem(U{1},Nq(2),1).*repmat(U{2},Nq(1),1) ; toc
norm(M1(:)-M2(:))


%% SPARSE MATRICES
clc ; clearvars ;
Nq = [10 20] ; 
fillRatio = 1e-2 ;
K = 10 ;

Nnz = round(fillRatio*Nq*K) ;
U = {} ; 
for mm=1:numel(Nq) 
    ii = randi(Nq(mm)*K,[Nnz(mm) 1]) ;
    vv = rand([Nnz(mm) 1]) ;
    U{mm} = reshape(sparse(ii,1,vv,Nq(mm)*K,1),[Nq(mm) K]) ;
end

tic ; M1 = pkg.math.khatrirao(U{:}) ; toc , M1
tic ; M2 = repelem(U{1},Nq(2),1).*repmat(U{2},Nq(1),1) ; toc , M2

end






