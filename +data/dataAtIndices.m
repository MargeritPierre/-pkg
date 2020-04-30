function VAL = dataAtIndices(DATA,IND)
% return values VAL corresponding to the data DATA queried and node indices IND
% IND [...szInd...] (any size, an index<0 is taken as non valid (==NaN))
% DATA [nNodes ...szData...]
% VAL [szInd szData]
    szInd = size(IND) ;
    szData = size(DATA) ;
    % Reshape the data
        nNodes = szData(1) ; 
        szData = szData(:,end) ;
        DATA = DATA(:,:) ;
    % Set the values
        valid = IND>0 & ~(IND>nNodes) ; % NaNs will return false also
        VAL = NaN([numel(IND) prod(szData)]) ;
        VAL(valid,:) = DATA(IND(valid),:) ; 
    % Reshape
        VAL = reshape(VAL,[szInd szData]) ;
end
