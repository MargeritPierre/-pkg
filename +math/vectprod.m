function V = vectprod(A,B,dim)
%VECTPROD return the vectorial product of two vectors
% A and B can be arrays; in this case, dim gives the dimension of the array
% along which to compute the product
% dim is by default the last common dimension of the two arrays
% constraint: size(A,dim) == size(B,dim) == {2,3}
% /!\ THE OUTPUT WILL ALWAYS BE 3D !!

    if nargin<2 ; error('Not enough arguments') ; end
    
    % Infos
    szA = size(A) ; szB = size(B) ;
    ndA = ndims(A) ; ndB = ndims(B) ;
    ndMin = min(ndA,ndB) ; ndMax = max(ndA,ndB) ;
    
    % Valid dimension: A singleton or B singleton or matching sizes
    validDim = szA(1:ndMin)==1 | szB(1:ndMin)==1 | szA(1:ndMin)==szB(1:ndMin) ;
    if any(~validDim) ; error('Wrong shape for A or B') ; end
    
    % Auto detect the dimension
    if nargin<3  
        validDim = szA(1:ndMin)<=3 & szB(1:ndMin)<=3 ;
        dim = find(validDim,1,'last') ;
    else
        if szA(dim)>3 || szB(dim)>3
            error('Cannot compute the vectorial product for more than 3 coordinates !') ;
        end
    end
    
    % Permute the input arrays to have dim as the first dimension
    permDims = [dim setdiff(1:ndMax,dim)] ;
    A = permute(A,permDims) ; 
    B = permute(B,permDims) ;
    
    % Pad arrays up to 3 coordinates
    A = padarray(A,[3-szA(dim) 0],0,'post') ;
    B = padarray(B,[3-szB(dim) 0],0,'post') ;
    
    % Compute the product
    V = A([2 3 1],:).*B([3 1 2],:) - A([3 1 2],:).*B([2 3 1],:) ;
    
    % Re-permute the output
    [~,permDims] = sort(permDims) ;
    V = permute(V,permDims) ;

end

