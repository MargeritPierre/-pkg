function [IN,iL] = findIn(values,list)
%FINDIN find values in a list
% values [nVal 1]
% list [nList 1]
% IN [nVal nList] sparse matrix
% memory-efficient implementation of IN = values(:)==list(:)' ;

nVal = numel(values) ;
nList = numel(list) ;

vv = [values(:)' list(:)'] ;
[uu,~,ic] = unique(vv) ; % vv=uu(ic) and uu = vv(ia)

isVal = sparse(1:nVal,ic(1:nVal),true,nVal,numel(uu)) ;
isList = sparse(1:nList,ic(nVal+1:end),true,nList,numel(uu)) ;

IN = logical(isVal*isList') ;

if nargout>1 ; [IN,iL] = find(IN) ; end

end


%% UNIT TESTS
function tests
%% 
values = 1:100000 ; list = 1:100 ;
values = randi(1e10,[1e5 1]) ; list = randi(1e10,[1e5 1]) ; 

disp(newline)
profile on
tic ; ON = values(:)==list(:)' ; toc
tic ; IN = pkg.data.findIn(values,list) ; toc
profile off

%sparse(ON)-IN




end