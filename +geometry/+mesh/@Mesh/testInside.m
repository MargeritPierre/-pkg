clc,clf,clearvars ;

N = 10 ;
x = linspace(0,1,N+1)' ;
mesh = pkg.geometry.mesh.Mesh(x) ;
mesh.extrude([0 1],N) ;
mesh.extrude([0 0 1],N) ;

if 0 % introduce random triangles
    toTri = logical(randi([0 1],[mesh.nElems 1])) ;
    mesh.Elems = [mesh.Elems.subpart(~toTri) ; mesh.Elems.subpart(toTri).simplex] ;
end

plot(mesh) ;
axis equal

%%

bbox = mesh.boundingBox ;
bbox = bbox + 0.25*[-1 ; 1].*range(bbox,1) ;

P = rand(10000,mesh.nCoord).*range(bbox,1) + min(bbox,[],1) ;
tol = 1e-9 ;

clf,clc ;

colors = jet(mesh.nElems) ;
colors = colors(randperm(mesh.nElems),:) ;

% pl = plot(mesh) ;
% pl.Faces.FaceColor = 'flat' ; 
% pl.Faces.FaceVertexCData = colors ;
% pl.Faces.FaceAlpha = 0.5 ;

axis equal
%plot(P(:,1),P(:,2),'.k','linewidth',.5) 


Xe = mesh.Elems.dataAtIndices(mesh.Nodes) ;

tic
p = permute(P,[1 3 2]) ;
Xmin = permute(min(Xe,[],2),[2 1 3])-tol ;
Xmax = permute(max(Xe,[],2),[2 1 3])+tol ;
in = all(Xmin<=p,3) & all(Xmax>=p,3) ;
toc

profile on
tic
bboxes = [min(Xe,[],2) max(Xe,[],2)] + [-1 1]*tol ;
[ip,ie] = pkg.data.inDomain(P,bboxes) ;
toc
profile viewer
%%
pOut = setdiff(1:size(P,1),ip) ;
scatter(P(pOut,1),P(pOut,2),20,'k','filled')
scatter(P(ip,1),P(ip,2),20,colors(ie,:),'filled')

