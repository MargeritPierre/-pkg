function list = sparse2list(M)
%SPARSE2LIST Convert a sparse (incidence) matrix to its connectivity list
%representation
% M: [nI nJ]
% list: [nI nMaxJJbyI]

% Find nonzero values
    [ii,jj,vv] = find(M) ;
    nRows = max(ii(:)) ;
    
% Sort by values (list column ordering)
    [~,ind] = sort(vv,'ascend') ;
    ii = ii(ind) ; jj = jj(ind) ;
    
% Sort row indices
    [ii,ind] = sort(ii,'ascend') ;
    jj = jj(ind) ;
    
% Number of valid indices by row
    ni = accumarray(ii(:),1,[nRows 1],[],0) ;
    
% Column indices
    cc = pkg.data.colon(ones(1,nRows),ni') ;
    
% Build the list
    list = zeros(nRows,max(ni)) ;
    list(sub2ind(size(list),ii(:),cc(:))) = jj(:) ;
    
end

