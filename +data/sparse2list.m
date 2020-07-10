function list = sparse2list(M)
%SPARSE2LIST Convert a sparse (incidence) matrix to its connectivity list
%reprezentation
% M: [nI nJ]
% list: [nI nMaxJJbyI]

    [ii,jj,vv] = find(M) ;
    nRows = max(ii(:)) ;
    
    [~,ind] = sort(vv,'ascend') ;
    ii = ii(ind) ; jj = jj(ind) ;
    
    elems = accumarray(ii(:),jj(:),[nRows 1],@(x){x},{0}) ;
    nEl = cellfun(@numel,elems) ;
    
    list = zeros(max(nEl),nRows) ;
    list(sub2ind(size(list),nEl(:)',1:nRows)) = 1 ;
    list = flip(cumsum(flip(list,1),1),1) ;
    list(logical(list)) = cat(1,elems{:})  ;
    
    list = list' ;
    
end

