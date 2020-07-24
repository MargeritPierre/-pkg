function d = det(varargin)
%DET Return the determinant of stacks of small matrices
% Applies element-wise operations if possible
% 
%   d = DET(A)
%       where size(A) = [n n sz1 sz2 sz3 ..]
%   d = DET(A11,A21,A31,...,An1,A12,...,Ann)


% Process inputs
if nargin==1 % INPUT : d = DET(A)
    A = varargin{1} ;
    % Deal with sizes
        SZ = size(A) ; 
        n = SZ(1) ; 
        if SZ(2)~=n ; error('Matrix must be square !') ; end
        SZ(1:2) = 1 ;
else % INPUT : d = DET(A11,A21,A31,...,An1,A12,...,Ann)
    n = sqrt(nargin) ; 
    % Deal with sizes
        if n~=round(n) ; error('Wrong number of arguments') ; end
        SZ = size(varargin{1}) ;
        A = cat(1,varargin{:}) ;
        A = reshape(A,[n n SZ]) ;
end

% COMPUTE THE DETERMINANT
switch n
    case 0
        d = ones([SZ 1 1]) ; 
    case 1
        d = A ; 
    case 2
        d = A(1,1,:).*A(2,2,:)-A(2,1,:).*A(1,2,:) ;
    case 3
        d = (A(2,1,:).*A(3,2,:)-A(3,1,:).*A(2,2,:)).*A(1,3,:) ...
            - (A(1,1,:).*A(3,2,:)-A(3,1,:).*A(1,2,:)).*A(2,3,:) ...
            + (A(1,1,:).*A(2,2,:)-A(2,1,:).*A(1,2,:)).*A(3,3,:) ...
            ;
    otherwise % n is too big: no analytical formula
        d = pkg.data.matrixfun(@det,A) ;
end
d = reshape(d,SZ) ;

return




%% TESTS

%% ONE INPUT d = det(A)
    s = 2 ; N = 10000 ;
    A = rand(s,s,N) ;
    tic ; d0 = pkg.data.matrixfun(@det,A) ; t0 = toc ;
    tic ; d = pkg.math.det(A) ; t = toc ;
    disp(['DET(A): ' ...
            's=' num2str(s) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(d(:)-d0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 
    
%% MULTIPLE INPUTS d = det(A)
    s = 3 ; N = 1 ;
    A = rand(s,s,N) ;
    tic ; d0 = pkg.data.matrixfun(@det,A) ; t0 = toc ;
    A = reshape(num2cell(A,3),[],1) ;
    tic ; d = pkg.math.det(A{:}) ; t = toc ;
    disp(['DET(A): ' ...
            's=' num2str(s) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(d(:)-d0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 

