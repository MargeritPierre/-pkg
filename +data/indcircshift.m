function B = indcircshift(A,K,dim)
%INDCIRCSHIFT circular shift of valid indices only


valid = A>0 ;

idx = cell(ndims(A),1) ; 
[idx{:}] = ind2sub(size(A),find(valid(:))) ;




end

