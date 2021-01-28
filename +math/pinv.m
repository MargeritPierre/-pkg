function varargout = pinv(varargin)
%PINV Return the pseudo-inverse of stacks of small matrices
% Applies element-wise operations if possible
% 
%   iA = PINV(A)
%       where size(A) = [n m sz1 sz2 sz3 ..]
%             size(iA) = [m n sz1 sz2 sz3 ..] !
%   iA = PINV([n m],A11,A21,A31,...,An1,A12,...,Anm)
%   [iA11,iA21,iA31,...,iAm1,iA12,...,iAmn] = PINV([n m],A11,A21,A31,...,An1,A12,...,Anm)


% Process inputs
if nargin==1 % INPUT : d = PINV(A)
    A = varargin{1} ;
    % Deal with sizes
        SZ = size(A) ; 
        n = SZ(1) ; m = SZ(2) ;
        SZ = [SZ(3:end) 1 1] ;
else % INPUT : d = INV([n m],A11,A21,A31,...,An1,A12,...,Ann)
    n = varargin{1}(1) ; m = varargin{1}(2) ;
    varargin = varargin(2:end) ;
    % Deal with sizes
        if nargin~=n*m+1 ; error('Wrong number of arguments') ; end
        SZ = size(varargin{1}) ;
        A = cat(1,varargin{:}) ;
        A = reshape(A,[n m SZ]) ;
end

% Is matrix square ?
if n==m 
    varargout = cell(1,nargout) ;
    [varargout{:}] = pkg.math.inv(A) ; 
    return ; 
end

% Check queried outputs
if nargout>1 && nargout~=n*m ; error('Wrong number of outputs') ; end

% COMPUTE THE PSEUDO-INVERSE
A = reshape(A,n,m,[]) ;
if m>n ; A = permute(A,[2 1 3]) ; end
switch m
    case 0
        iA = [] ; 
    case 1 % pinv(v) = conj(v')/norm(v')^2
        iA = permute(conj(A),[2 1 3])./sum(abs(A).^2,1) ; 
    otherwise % pinv(A) = inv(A'*A)*A'
        At = permute(conj(A),[2 1 3]) ;
        AxA = pkg.math.mtimes(At,A) ;
        iAxA = pkg.math.inv(AxA) ;
        iA = pkg.math.mtimes(iAxA,At) ;
end
if m>n ; iA = permute(iA,[2 1 3]) ; end

% Process for output
iA = reshape(iA,[m n SZ]) ;
if nargout<=1
    varargout = {iA} ;
else
    iA = permute(iA,[2+(1:numel(SZ)) 1 2]) ;
    varargout = num2cell(iA,1:numel(SZ)) ;
end

return


%% TESTS

%% ONE INPUT iA = pinv(A)
    n = 3 ; m = 3 ; N = 100000 ;
    A = rand(n,m,N)  + 1i*rand(n,m,N) ;
    tic ; iA0 = pkg.data.matrixfun(@pinv,A) ; t0 = toc ;
    tic ; iA = pkg.math.pinv(A) ; t = toc 
    disp(['PINV(A): ' ...
            's=[' num2str(n) ',' num2str(m) ']' ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(iA(:)-iA0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 
    
%% MULTIPLE INPUTS iA = inv(A11,A21,...,Anm)
    n = 3 ; m = 2 ; N = 10000 ;
    A = rand(n,m,N)  + 1i*rand(n,m,N) ;
    tic ; iA0 = pkg.data.matrixfun(@pinv,A) ; t0 = toc ;
    A = reshape(num2cell(A,3),[],1) ;
    tic ; iA = pkg.math.pinv([n m],A{:}) ; t = toc ;
    disp(['INV(A): ' ...
            's=' num2str(n) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(iA(:)-iA0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 
    
%% MULTIPLE OUTPUTS [iA{:}] = inv(A)
    n = 3 ; m = 2 ; N = 10000 ;
    A = rand(n,m,N)  + 1i*rand(n,m,N) ;
    tic ; iA0 = pkg.data.matrixfun(@pinv,A) ; t0 = toc ;
    iA = cell(n*m,1) ;
    tic ; [iA{:}] = pkg.math.pinv(A) ; t = toc ;
    iAr = reshape(cat(1,iA{:}),[n m N]) ;
    disp(['INV(A): ' ...
            's=' num2str(n) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(iAr(:)-iA0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 
    
%% MULTIPLE INPUTS && OUTPUTS [iA{:}] = inv(A11,A21,...,Ann)
    n = 3 ; m = 2 ; N = 10000 ;
    A = rand(n,m,N)  + 1i*rand(n,m,N) ;
    tic ; iA0 = pkg.data.matrixfun(@pinv,A) ; t0 = toc ;
    iA = cell(n*m,1) ;
    A = reshape(num2cell(A,3),[],1) ;
    tic ; [iA{:}] = pkg.math.pinv([n m],A{:}) ; t = toc ;
    iAr = reshape(cat(1,iA{:}),[n m N]) ;
    disp(['INV(A): ' ...
            's=' num2str(n) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(iAr(:)-iA0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 

