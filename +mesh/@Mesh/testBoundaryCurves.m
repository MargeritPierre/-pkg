clearvars -except mesh

edg = mesh.Edges(mesh.boundaryEdges,:) ;
edg = sort(edg,2) ;
[ee,ii] = sort(edg(:),1) ;
[eee,nnn] = ind2sub(size(edg),ii) ;
isflipped = nnn==2 ;
crv = edg(eee,:) ;

[crv,ia,ic] = unique(crv,'rows')

%%
% [~,M_IDX] = ismember(E(:,1),E(:,2));
% [~,S_IDX] = sort(M_IDX);
% G = E(S_IDX,:)

clc
clearvars -except mesh

edg = mesh.Edges(mesh.boundaryEdges,:) ;
edg = sort(edg,2) ;

% Switch edges if needed
% check that an indice is not present two times in a column
for col = 1:2
    [~,ind] = ismember(edg(:,col),edg(:,col)) ; % find index of the first equal value
    toSwitch = ind(:)'~=1:numel(ind) ; % if it is not the same index, switch
    edg(toSwitch,:) = flip(edg(toSwitch,:),2) ;
end

%% Find the next edge
[~,ind] = ismember(edg(:,2),edg(:,1)) ;
%[~,ind] = sort(ind) ;
%
%crv1 = edg(ind1(valid1,:),:)

%crv = edg*NaN ;
%crv(

%%
clc
clearvars -except mesh

edg = mesh.Edges(mesh.boundaryEdges,:) ;

[nod,ia,ic] = unique(edg(:)) ;
nod = 1:numel(nod) ;
edg = reshape(nod(ic),[],2) ;

g = graph(edg(:,1),edg(:,2)) ;
crv = conncomp(g,'OutputForm','cell') ;

%crv(end+1,:) = {NaN} ;
crv = cat(2,crv{:}) ;
crv = crv(1:end-1) ;

clf
plot(mesh.Nodes(crv,1),mesh.Nodes(crv,2))

%%
clc
clearvars -except mesh
tic

edg = mesh.Edges(mesh.boundaryEdges,:) ;
e2n = mesh.edge2nod ;
n2n = e2n*e2n' ;
[n1,n2] = find(n2n==1) ;


toc






