function [k,U] = sort(k,U,dim,realPos,gammaMax) 
% Sort wavenumbers & eigenvectors
if nargin<3 || iempty(dim) ; dim = 1 ; end
dimU = dim+1 ;

[k,isort] = sort(k,dim,'ascend','ComparisonMethod','real') ;

isort = permute(isort,[dim 1:dim-1 dim+1:ndims(isort)]) ;
permU = [dimU 1:dimU-1 dimU+1:ndims(U)] ;
U = permute(U,permU) ;

for ii = 1:numel(isort)/size(isort,1) 
    U(:,:,ii) = U(isort(:,ii),:,ii) ;
end

[~,iperm] = sort(permU) ;
U = permute(U,iperm) ;


end