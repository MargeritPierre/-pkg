function VAL = meanDataAtIndices(DATA,IND,DIM)
% Return the mean of DATA evaluated at indices IND in the dimension DIM
    if nargin<3 ; DIM = ndims(IND) ; end
    VAL = mean(pkg.data.dataAtIndices(DATA,IND),DIM,'omitnan') ;
end