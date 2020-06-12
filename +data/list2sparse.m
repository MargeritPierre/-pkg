function M = list2sparse(list,values)
% LIST2SPARSE Convert a connectivity list to its sparse representation
% Corresponds to the incidence matrix
% Use "values" input argument to put custom values in the matrix

    if nargin<2 ; values = 'indices' ; end
    
    if isnumeric(values) % User customized values
        vvv = values + zeros(size(list)) ;
    else
        switch values
            case 'logical' % M(ii,jj) = true if ismember(ii,list(jj,:))
                vvv = true(size(list)) ;
            case 'indices' % M(ii,jj) = find(list(jj,:)==ii) ;
                vvv = repmat(1:size(list,2),[size(list,1) 1]) ;
            case 'mean' % Used to compute the mean of "nodal" values on each "element"
                vvv =  true(size(list)) ;
        end
    end
    
    % Integers are the indice of the node in the elem list
    valid = list>0 ;
    eee = repmat((1:size(list,1))',[1 size(list,2)]) ;
    M = sparse(double(list(valid)),eee(valid),vvv(valid)) ;
    
    % Eventually compute the mean
    if strcmp(values,'mean') ; M = M.*(1./sum(M,2)) ; end
    
end

