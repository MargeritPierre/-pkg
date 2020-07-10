function [C,ind] = combinations(varargin)
%COMBINATIONS return all possible combinations of a set of possible values
% varargin = {[v11 v12 v13 ..] [v21 v22 ..] ..}
% C = [nComb nInputs] values
% ind = [nComb nInputs] corresponding indices in each bin

    if numel(varargin)==1 && iscell(varargin{1})
        varargin = varargin{1} ;
    end
    
    varargin = reshape(varargin,1,[]) ;

    nInputs = numel(varargin) ;
    nBinsByInput = cellfun(@numel,varargin) ;
    
    if any(nBinsByInput==0) ; C = [] ; ind = [] ; return ; end
    
    ind = arrayfun(@colon,nBinsByInput*0+1,nBinsByInput,'UniformOutput',false) ;
    [ind{:}] = ndgrid(ind{:}) ;
    ind = cat(nInputs+1,ind{:}) ;
    ind = reshape(ind,[],nInputs) ;
    
    C = cat(2,varargin{:}) ;
    linInd = ind + cumsum([0 nBinsByInput(1:end-1)]) ;
    C = C(linInd) ;
end

