function C = innerprod(A,B,varargin)
%INNERPROD return the inner product of two "tensors" A and B (arrays...)
% varargin: either one ore two OPTIONNAL arguments:
% 	- one argument: ORDER (opt., default: 1) allows to compute double (2), triple (3), etc. inner products 
%       in this case the summation is done along the N last dimensions of A
%       and first dimensions of B
%   - two arguments: dimsA and dimsB: 
%       specifies on which dimensions the summation is done
%       dimsA and dimsB must contain the same number of elements

    if nargin<2 ; error('Not enough arguments') ; end
    
    % Infos
    szA = size(A) ;
    szB = size(B) ;
    ndA = numel(szA) ;
    ndB = numel(szB) ;
    
    % Input process
    switch numel(varargin)
        case 0
            dimsA = ndA ;
            dimsB = ndB ;
        case 1
            order = varargin{1} ;
            dimsA = ndA:-1:ndA-order+1 ;
            dimsB = 1:order ;
        case 2
            dimsA = varargin{1} ;
            dimsB = varargin{2} ;
            if numel(dimsA)~=numel(dimsB)
                error('Number of contracting dimensions must match') ;
            end
        otherwise
            error('Too many input arguments') ;
    end
    
    % Check dimensions
    validDims = szA(dimsA)==szB(dimsB) ;
    if any(~validDims)
        error('Wrong shape for A and B') ;
    end
    
    % Permute the input arrays
    otherDimsA = setdiff(1:ndA,dimsA) ;
    otherDimsB = setdiff(1:ndB,dimsB) ;
    A = permute(A,[otherDimsA dimsA]) ;
    B = permute(B,[dimsB otherDimsB]) ;
    
    % To matrix (for product vectorization)
    A = reshape(A,[prod(szA(otherDimsA)) prod(szA(dimsA))]) ;
    B = reshape(B,[prod(szB(dimsB)) prod(szB(otherDimsB))]) ;
    
    % Compute the product using the optimized built-in method
    C = A*B ;
    
    % Reshape the output
    C = reshape(C,[szA(otherDimsA) szB(otherDimsB) 1 1]) ;

end

