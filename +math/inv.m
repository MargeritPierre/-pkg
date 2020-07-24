function varargout = inv(varargin)
%INV Return the inverse of stacks of small matrices
% Applies element-wise operations if possible
% 
%   iA = INV(A)
%       where size(A) = [n n sz1 sz2 sz3 ..]
%   iA = INV(A11,A21,A31,...,An1,A12,...,Ann)
%   [iA11,iA21,iA31,...,iAn1,iA12,...,iAnn] = INV(A11,A21,A31,...,An1,A12,...,Ann)


% Process inputs
if nargin==1 % INPUT : d = INV(A)
    A = varargin{1} ;
    % Deal with sizes
        SZ = size(A) ; 
        n = SZ(1) ; 
        if SZ(2)~=n ; error('Matrix must be square !') ; end
        SZ(1:2) = [] ;
else % INPUT : d = INV(A11,A21,A31,...,An1,A12,...,Ann)
    n = sqrt(nargin) ; 
    % Deal with sizes
        if n~=round(n) ; error('Wrong number of arguments') ; end
        SZ = size(varargin{1}) ;
        A = cat(1,varargin{:}) ;
        A = reshape(A,[n n SZ]) ;
end

% Check queried outputs
    if nargout>1 && nargout~=n^2 ; error('Wrong number of outputs') ; end

% COMPUTE THE INVERSE
switch n
    case 0
        iA = [] ; 
    case 1
        iA = 1./A ; 
    case 2
        d = A(1,1,:).*A(2,2,:)-A(2,1,:).*A(1,2,:) ;
        iA = (1./d).*[A(2,2,:) -A(1,2,:) ; -A(2,1,:) A(1,1,:)] ;
    case 3
        % Determinant..
        d = (A(2,1,:).*A(3,2,:)-A(3,1,:).*A(2,2,:)).*A(1,3,:) ...
            - (A(1,1,:).*A(3,2,:)-A(3,1,:).*A(1,2,:)).*A(2,3,:) ...
            + (A(1,1,:).*A(2,2,:)-A(2,1,:).*A(1,2,:)).*A(3,3,:) ...
            ;
        % Cofactors ..
        iA = zeros([n n prod(SZ)]) ;
        iA(1,1,:) = A(2,2,:).*A(3,3,:)-A(2,3,:).*A(3,2,:) ;
        iA(1,2,:) = A(2,1,:).*A(3,3,:)-A(2,3,:).*A(3,1,:) ;
        iA(1,3,:) = A(2,1,:).*A(3,2,:)-A(2,2,:).*A(3,1,:) ;
        iA(2,1,:) = A(1,2,:).*A(3,3,:)-A(1,3,:).*A(3,2,:) ;
        iA(2,2,:) = A(1,1,:).*A(3,3,:)-A(1,3,:).*A(3,1,:) ;
        iA(2,3,:) = A(1,1,:).*A(3,2,:)-A(1,2,:).*A(3,1,:) ;
        iA(3,1,:) = A(1,2,:).*A(2,3,:)-A(1,3,:).*A(2,2,:) ;
        iA(3,2,:) = A(1,1,:).*A(2,3,:)-A(1,3,:).*A(2,1,:) ;
        iA(3,3,:) = A(1,1,:).*A(2,2,:)-A(1,2,:).*A(2,1,:) ;
        iA = (1./d(:,:,:)).*(permute(iA,[2 1 3]).*[1 -1 1 ; -1 1 -1 ; 1 -1 1]) ; 
    otherwise % n is too big: no analytical formula
        iA = pkg.data.matrixfun(@inv,A) ;
end

% Process for output
iA = reshape(iA,[n n SZ]) ;
if nargout<=1
    varargout = {iA} ;
else
    iA = permute(iA,[2+(1:numel(SZ)) 1 2]) ;
    varargout = num2cell(iA,1:numel(SZ)) ;
end

return



%% OPTIMIZED ACCESS TO MEMORY: 

% A = permute(A(:,:,:),[3 2 1]) ;
% % Determinant..
% d = (A(:,2,1).*A(:,3,2)-A(:,3,1).*A(:,2,2)).*A(:,1,3) ...
%     - (A(:,1,1).*A(:,3,2)-A(:,3,1).*A(:,1,2)).*A(:,2,3) ...
%     + (A(:,1,1).*A(:,2,2)-A(:,2,1).*A(:,1,2)).*A(:,3,3) ...
%     ;
% % Cofactors ..
% in = zeros([prod(SZ) n n]) ;
% in(:,1,1) = A(:,2,2).*A(:,3,3)-A(:,2,3).*A(:,3,2) ;
% in(:,1,2) = A(:,2,1).*A(:,3,3)-A(:,2,3).*A(:,3,1) ;
% in(:,1,3) = A(:,2,1).*A(:,3,2)-A(:,2,2).*A(:,3,1) ;
% in(:,2,1) = A(:,1,2).*A(:,3,3)-A(:,1,3).*A(:,3,2) ;
% in(:,2,2) = A(:,1,1).*A(:,3,3)-A(:,1,3).*A(:,3,1) ;
% in(:,2,3) = A(:,1,1).*A(:,3,2)-A(:,1,2).*A(:,3,1) ;
% in(:,3,1) = A(:,1,2).*A(:,2,3)-A(:,1,3).*A(:,2,2) ;
% in(:,3,2) = A(:,1,1).*A(:,2,3)-A(:,1,3).*A(:,2,1) ;
% in(:,3,3) = A(:,1,1).*A(:,2,2)-A(:,1,2).*A(:,2,1) ;
% in = (1./d).*in ;
% in = permute(in,[2 3 1]).*[1 -1 1 ; -1 1 -1 ; 1 -1 1] ;




%% TESTS

%% ONE INPUT iA = inv(A)
    s = 3 ; N = 10000 ;
    A = rand(s,s,N)  + 1i*rand(s,s,N) ;
    tic ; iA0 = pkg.data.matrixfun(@inv,A) ; t0 = toc ;
    tic ; iA = pkg.math.inv(A) ; t = toc ;
    disp(['INV(A): ' ...
            's=' num2str(s) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(iA(:)-iA0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 
    
%% MULTIPLE INPUTS iA = inv(A11,A21,...,Ann)
    s = 3 ; N = 100000 ;
    A = rand(s,s,N) ;
    tic ; iA0 = pkg.data.matrixfun(@inv,A) ; t0 = toc ;
    A = reshape(num2cell(A,3),[],1) ;
    tic ; iA = pkg.math.inv(A{:}) ; t = toc ;
    disp(['INV(A): ' ...
            's=' num2str(s) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(iA(:)-iA0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 
    
%% MULTIPLE OUTPUTS [iA{:}] = inv(A)
    s = 3 ; N = 1000 ;
    A = rand(s,s,N) ;
    tic ; iA0 = pkg.data.matrixfun(@inv,A) ; t0 = toc ;
    iA = cell(s^2,1) ;
    tic ; [iA{:}] = pkg.math.inv(A) ; t = toc ;
    iAr = reshape(cat(1,iA{:}),[s s N]) ;
    disp(['INV(A): ' ...
            's=' num2str(s) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(iAr(:)-iA0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 
    
%% MULTIPLE INPUTS && OUTPUTS [iA{:}] = inv(A11,A21,...,Ann)
    s = 3 ; N = 1000 ;
    A = rand(s,s,N) ;
    tic ; iA0 = pkg.data.matrixfun(@inv,A) ; t0 = toc ;
    iA = cell(s^2,1) ;
    A = reshape(num2cell(A,3),[],1) ;
    tic ; [iA{:}] = pkg.math.inv(A{:}) ; t = toc ;
    iAr = reshape(cat(1,iA{:}),[s s N]) ;
    disp(['INV(A): ' ...
            's=' num2str(s) ...
            ', N=' num2str(N) ...
            ', Error: ' num2str(norm(iAr(:)-iA0(:))) ...
            ', speed ratio: ' num2str(t/t0) ...
        ]) 

