function varargout = matrixfun(fun,varargin)
%MATRIXFUN Apply a function on a stack of matrices
% fun is a function handle
% varargin{i} contains the array Ai of size [nRows_i nColumns_i SZ]
% SZ = [sz1 sz2 ...] is [1] or a constant over the input matrices
% UniformOutput is true by default

% UniformOutput being set ?
uoArg = strcmp(varargin,'UniformOutput') ;
if any(uoArg) % 'UniformOuput' has been provided
    uoArg = find(uoArg) ;
    UniformOutput = varargin(uoArg+1) ;
    varargin(uoArg+(0:1)) = [] ;
else
    UniformOutput = true ;
end

% Deal with sizes..
    SZ = cellfun(@size,varargin,'UniformOutput',false) ;
    NDIMS = cellfun(@numel,SZ) ;
% Take the biggest array as reference (not general...)
    [~,biggest] = max(NDIMS) ;
    SZ = SZ{biggest}(3:end) ;
    N = prod(SZ) ;
    
% Init input & output
    out = cell(nargout,N) ;
    
% PROCESS
    if nargin-1==1 && nargout==1 
    % one in & one out
        in = varargin{1} ;
        for nn = 1:N
            out{nn} = fun(in(:,:,nn)) ;
        end
    elseif nargin-1==1 
    % one in & multiple out
        in = varargin{1} ;
        for nn = 1:N
            [out{:,nn}] = fun(in(:,:,nn)) ;
        end
    elseif nargout==1 
    % multiple in & one out
        in = cell(nargin-1,1) ;
        for nn = 1:N
            for ii = 1:nargin-1 ; in{ii} = varargin{ii}(:,:,nn) ; end
            out{nn} = fun(in{:}) ;
        end
    else
    % multiple in & multiple out
        in = cell(nargin-1,1) ;
        for nn = 1:N
            for ii = 1:nargin-1 ; in{ii} = varargin{ii}(:,:,nn) ; end
            [out{:,nn}] = fun(in{:}) ;
        end
    end

% Apply to each matrix of the array

% Uniform output ?
    varargout = cell(nargout,1) ;
    if UniformOutput
        for oo = 1:nargout
            varargout{oo} = reshape(cat(3,out{oo,:}),[size(out{oo,1}) SZ 1 1]) ;
        end
    else
        for oo = 1:nargout
            varargout{oo} = reshape(varargout{oo,:},[SZ 1 1]) ;
        end
    end

return



%% TESTS

%% DETERMINANT d = det(A)
    s = 100 ; N = 1000 ;
    A = rand(s,s,N) ;
    tic ; d0 = zeros(1,1,N) ; for nn = 1:N ; d0(nn) = det(A(:,:,nn)) ; end ; t0 = toc ;
    tic ; d = pkg.data.matrixfun(@det,A) ; t = toc ;
    disp(['DET: ' ...
            's=' num2str(s) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(d(:)-d0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 
    
%% INVERSE i = inv(A)
    s = 3 ; N = 1000 ;
    A = rand(s,s,N) ;
    tic ; i0 = zeros(s,s,N) ; for nn = 1:N ; i0(:,:,nn) = inv(A(:,:,nn)) ; end ; t0 = toc ;
    tic ; i = pkg.data.matrixfun(@inv,A) ; t = toc ;
    disp(['INV: ' ...
            's=' num2str(s) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(i(:)-i0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 
    
%% MLDIVIDE A.x = b
    s = 2 ; N = 10000 ;
    A = rand(s,s,N) ; b = rand(s,1,N) ;
    tic ; 
        x0 = zeros(s,1,N) ; d0 = zeros(s,s,N) ; 
        %for nn = 1:N ; x0(:,:,nn) = mldivide(A(:,:,nn),b(:,:,nn)) ; end 
        for nn = 1:N ; x0(:,:,nn) = A(:,:,nn)\b(:,:,nn) ; end 
    t0 = toc ;
    tic ; x = pkg.data.matrixfun(@mldivide,A,b) ; t = toc ;
    disp(['MLDIVIDE: ' ...
            's=' num2str(s) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(x(:)-x0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 
    
%% Generalized EigenValues [W,D] = eig(A,B)
    s = 20 ; N = 100 ;
    A = rand(s,s,N) ; B = rand(s,s,N) ;
    tic ; 
        w0 = zeros(s,s,N) ; d0 = zeros(s,s,N) ; 
        for nn = 1:N ; [w0(:,:,nn),d0(:,:,nn)] = eig(A(:,:,nn),B(:,:,nn)) ; end ; 
    t0 = toc ;
    tic ; [w,d] = pkg.data.matrixfun(@eig,A,B) ; t = toc ;
    disp(['EIG: ' ...
            's=' num2str(s) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(w(:)-w0(:)) + norm(d(:)-d0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 







