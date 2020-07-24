clc,clf,clearvars ;

N = 3 ;
x = linspace(0,1,N+1)' ;
mesh = pkg.geometry.mesh.Mesh(x) ;
mesh.extrude([0 1],N) ;
%mesh.extrude([0 0 1],N) ;

if 0 % introduce random triangles
    toTri = logical(randi([0 1],[mesh.nElems 1])) ;
    mesh.Elems = [mesh.Elems.subpart(~toTri) ; mesh.Elems.subpart(toTri).simplex] ;
end
if 1 % distort the grid a bit..
    mesh.Nodes = mesh.Nodes + 0*rand(mesh.nNodes,mesh.nCoord)*(1/N) ;
end

plot(mesh) ;
axis equal

%% TEST IF POINTS ARE INSIDE ELEMENTS

bbox = mesh.boundingBox ;
bbox = bbox + 0.05*[-1 ; 1].*range(bbox,1) ;
P = rand(1000000,mesh.nCoord).*range(bbox,1) + min(bbox,[],1) ;
tol = 1e-9 ;

clc,clf ;
pl = plot(mesh) ;
%plot(P(:,1),P(:,2),'.k','linewidth',.5) 

profile on
tic
[M,ip,ie] = mesh.isInside(P,mesh.Elems) ;
toc
profile viewer

colors = jet(mesh.nElems) ;
colors = colors(randperm(mesh.nElems),:) ;

pl.Faces.FaceColor = 'flat' ; 
pl.Faces.FaceVertexCData = colors ;
pl.Faces.FaceAlpha = 0.5 ;

axis equal

pOut = setdiff(1:size(P,1),ip) ;
scatter(P(pOut,1),P(pOut,2),20,'k','filled')
scatter(P(ip,1),P(ip,2),20,colors(ie,:),'filled')

%%
clc
tic
M = mesh.interpMat(P) ;
toc


%%
clc
nCoord = 3 ;
nP = 100 ;
nDom = 100 ;
time = false ;

P = rand(nP,nCoord) ;
domains = sort(rand(nDom,2,nCoord),2) ; % Random domains: lot of overlapping !
domains = permute(repmat(reshape(repelem(linspace(0,1,nDom+1),[1 2*ones(1,nDom-1) 1]),[2 nDom]),[1 1 nCoord]),[2 1 3]) ; % well-separated domains: no overlap

tic
p = permute(P,[1 3 2]) ;
dom = permute(domains,[2 1 3]) ;
in = all(dom(1,:,:)<=p,3) & all(dom(2,:,:)>=p,3) ;
%[ip,ie] = find(in) ;
toc

if time ; profile on ; end
tic
[ip,ie] = pkg.data.inDomain(P,domains) ;
toc
if time ; profile viewer ; end
