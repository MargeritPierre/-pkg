function list = sparse2list(M)
%SPARSE2LIST Convert a sparse (incidence) matrix to its connectivity list
%reprezentation

    [ii,jj,vv] = find(M) ;
    nRows = max(jj(:)) ;
    nCols = sum(logical(M),1) ;
    
    [~,ind] = sort(vv,'ascend') ;
    ii = ii(ind) ; jj = jj(ind) ; vv = vv(ind) ;
    
    elems = accumarray(jj(:),ii(:),[nRows 1],@(x){x},{0}) ;
    nEl = cellfun(@numel,elems) ;
    
    list = zeros(nRows,max(nEl)) ;
    list(sub2ind(size(list),1:nRows,nEl(:)')) = 1 ;
    list = flip(cumsum(flip(list,2),2),2) ;
    list(logical(list)) = cat(1,elems{:})  ;
    
end

