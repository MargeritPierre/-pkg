%clearvars
cla ; axis equal

% SHAPE CONTOUR DEFINITION
Nc = 1001 ;
if 0 % rectangle
    Cp = [-1 -1 ; 1 -1 ; 1 1 ; -1 1].*[10 8] ;
    Cp = interp1(Cp([1:end,1],:),linspace(1,size(Cp,1)+1,Nc)) ;
    Cp = Cp(1:end-1,:) ;
elseif 0 % ellipse
    a = 10 ; b = 5 ;
    phi = -18*pi/90 ;
    theta = linspace(0,2*pi,Nc) ;
    Cp = [a*cos(theta(:)) b*sin(theta(:))]*[cos(phi) -sin(phi) ; sin(phi) cos(phi)] ;
    Cp = Cp(1:end-1,:) ;
elseif 0 % rectangle with elliptical hole (not working)
    rect = [-1 -1 ; 1 -1 ; 1 1 ; -1 1].*[10 8] ;
    rect = interp1(rect([1:end,1],:),linspace(1,size(rect,1)+1,Nc)) ;
    a = 5 ; b = 5 ; phi = -18*pi/90 ; theta = linspace(0,2*pi,Nc) ;
    elli = [a*cos(theta(:)) b*sin(theta(:))]*[cos(phi) -sin(phi) ; sin(phi) cos(phi)] ;
    Cp = [rect ; elli] ;
elseif 1 % smooth shapee
    A = [16 0 6 6 0 0 0 0 0 0 0 0 0 2 0] ;
    phi = [0 0 0 0 0 -18 0 -3 37 0 0 0 0 -5 4]*pi/30 + 9*pi/90 ;
    theta = linspace(0,2*pi,Nc) ;
    R = sum(A(:).*cos(theta.*(0:numel(A)-1)' + phi(:)),1) ;
    Cp = R(:).*[cos(theta(:)) sin(theta(:))] ;
    Cp = Cp(1:end-1,:) ;
end

skel = pkg.geometry.Skeleton(Cp) ;
plot(Cp([1:end,1],1),Cp([1:end,1],2),'k','linewidth',1)
plot(skel)

% if 1 % skeleton
%     skel = pkg.geometry.Skeleton(Cp) ;
%     patches = skel.quadPatches ;
%     plot(skel.DelaunayMesh,'ShowLabels','Elems')%,'VisibleEdges','none') ;
%     plot(skel,'Color','r','EdgeWidth',1,'VisibleNodes','all','ShowLabels','Nodes')%,'HighlightEndNodes',true) ;
%     %plot(patches,'Color','b','EdgeWidth',0.5) ;
% elseif 1 % Simplex mesh
%     mesh = pkg.geometry.mesh.fromBoundary(Cp) ;
%     pl = plot(mesh) ;
% end

%%
mesh.LaplacianSmoothing([1 0]*0.8,100)
%mesh.Elems.Indices = padarray(delaunay(mesh.Nodes),[0 1],1,'pre') ;
%mesh.Elems = mesh.Elems.subpart(isInside(mesh.centroid)) ;
pl.update


%%

% DELAUNAY MESH TRIANGULATION
tri = delaunay(Cp) ;
mesh = pkg.geometry.mesh.Mesh('Nodes',Cp,'Elems',tri) ;
cen = mesh.centroid ;
mesh.removeElems(~inpolygon(cen(:,1),cen(:,2),Cp(:,1),Cp(:,2))) ; % Delete polygons outside
plot(mesh,'VisibleFaces','none')



% Elem circumcenters
    circum = pkg.math.circumcenter(permute(mesh.Elems.dataAtIndices(mesh.Nodes),[2 3 1])) ;
% Inner edges, that will be used to compute middle skeleton
    innerEdges = ~mesh.boundaryEdges ;
% Mesh corners
    corners = mesh.endNodes ;
% Center triangles, linked to only inner edges
    centerTris = (innerEdges'*mesh.elem2edge)==3 ;
% Circumcenter of center triangles
    S = circum(centerTris,:) ;
% Inner edge 2-face middle point
    face2edge = mesh.face2edge ;
    meanOverFaces = face2edge.*(1./sum(face2edge,2)) ;
    P = meanOverFaces(innerEdges,:)*circum ;
% Mesh corner nodes 
    C = mesh.Nodes(corners,:) ;
plot(circum(:,1),circum(:,2),'.m','markersize',10) ;
plot(S(:,1),S(:,2),'.r','markersize',20) ;
plot(C(:,1),C(:,2),'.b','markersize',20) ;


%% DEBUG MODE
cla ;
% Plot the mesh
plot(mesh,'VisibleFaces','none','VisibleNodes','all','ShowLabels',{'Nodes','Elems'}) ;
%plot(pkg.geometry.mesh.Mesh('Nodes',circum),'VisibleNodes','all','ShowLabels','Nodes')

