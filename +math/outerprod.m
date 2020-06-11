function C = outerprod(A,B,dimA,dimB)
%OUTERPROD return the outer product of two "tensors" A and B (arrays...)
% /!\ default behavior (without inputs dimA & dimB) 
%       - last singleton dimensions of A are not taken into account !
%       - first singleton dimensions of B are not taken into account !

    % Infos
    szA = size(A) ;
    szB = size(B) ;
    
    % Process inputs
    if nargin<2 ; error('Not enough arguments') ; end
    if nargin<3
        dimA = find(szA>1,1,'last') ; % last singleton dimensions of A
    end
    if nargin<4 
        dimB = find(szB>1,1,'first') ;  % first singleton dimensions of B
    end
    
    % Cut sizes
    szA = szA(1:dimA) ;
    szB = szB(dimB:end) ;
    
    % Reshape B
    B = reshape(B,[ones(1,numel(szA)) szB]) ;
    
    % Hadamard product
    C = A.*B ;
    
end

