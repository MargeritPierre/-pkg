function VAL = meanDataAtIndices(DATA,IND,DIM)
% Return the mean of DATA evaluated at indices IND in the dimension DIM
    if nargin<3 ; DIM = ndims(IND) ; end
    VAL = mean(pkg.data.dataAtIndices(DATA,IND),DIM,'omitnan') ;
    sz = size(VAL) ;
    VAL = reshape(VAL,[sz(setdiff(1:ndims(VAL),DIM)) 1]) ;
end