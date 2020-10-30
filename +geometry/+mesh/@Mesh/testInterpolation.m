clc,clf,clearvars ;

N = 3 ;
L = 0.1 ;
x = linspace(0,1,N+1)'*L ;
mesh = pkg.geometry.mesh.Mesh(x) ;
mesh.extrude([0 1]*L,N) ;
mesh.extrude([0 0 1]*L,N) ;

if 0 % introduce random triangles
    toTri = logical(randi([0 1],[mesh.nElems 1])) ;
    mesh.Elems = [mesh.Elems.subpart(~toTri) ; mesh.Elems.subpart(toTri).simplex] ;
end
if 1 % distort the grid a bit..
    mesh.Nodes = mesh.Nodes + .25*rand(mesh.nNodes,mesh.nCoord)*(L/N) ;
end

plot(mesh) ;
axis equal


%% COMPUTE THE GRADIENT OF A FUNCTION
x = mesh.Nodes ; % the space variable
tic ; G = mesh.gradMat ; toc % gradient matrices s.t. df/dx_i = G{i}*f ;
% Tests {{f} {df_dx1} {df_dx2}]}
tests = {} ;
tests(end+1,:) = {x(:,1) ones(mesh.nElems,1) zeros(mesh.nElems,1) zeros(mesh.nElems,1)} ;
tests(end+1,:) = {x(:,2) zeros(mesh.nElems,1) ones(mesh.nElems,1) zeros(mesh.nElems,1)} ;
tests(end+1,:) = {2*x(:,1)-5*x(:,2) + 1i*x(:,3) 2*ones(mesh.nElems,1) -5*ones(mesh.nElems,1) 1i*ones(mesh.nElems,1)} ;
tests(end+1,:) = {sum(x.^2,2) 2*mesh.Elems.meanDataAtIndices(x(:,1)) 2*mesh.Elems.meanDataAtIndices(x(:,2)) 0*mesh.Elems.meanDataAtIndices(x(:,2))} ;
% Apply
errors = NaN(size(tests,1),mesh.nCoord) ;
for tt = 1:size(tests,1)
    f = tests{tt,1} ;
    for cc = 1:min(mesh.nCoord,size(tests,2)-1)
        df_dxc = G{cc}*f ;
        errors(tt,cc) = norm(tests{tt,cc+1}-df_dxc) ;
    end
end
errors


%% COMPUTE THE SECOND GRADIENT OF A FUNCTION
x = mesh.Nodes ; % the space variable
tic ; G2 = mesh.grad2Mat ; toc % gradient matrices s.t. df/(dx_i.dx_j) = G2{i,j}*f ;
% Tests {{f} {df_dx1dx1} {df_dx1dx2} {df_dx2dx1} {df_dx2dx2}]}
tests = {} ;
tests(end+1,:) = {x(:,1) zeros(mesh.nElems,1) zeros(mesh.nElems,1) zeros(mesh.nElems,1) zeros(mesh.nElems,1)} ;
tests(end+1,:) = {x(:,2) zeros(mesh.nElems,1) zeros(mesh.nElems,1) zeros(mesh.nElems,1) zeros(mesh.nElems,1)} ;
tests(end+1,:) = {x(:,1).^2 2*ones(mesh.nElems,1) zeros(mesh.nElems,1) zeros(mesh.nElems,1) zeros(mesh.nElems,1)} ;
tests(end+1,:) = {x(:,1).*x(:,2) zeros(mesh.nElems,1) ones(mesh.nElems,1) ones(mesh.nElems,1) zeros(mesh.nElems,1)} ;
tests(end+1,:) = {x(:,1).*x(:,2).^2 zeros(mesh.nElems,1) 2*mesh.Elems.meanDataAtIndices(x(:,2)) 2*mesh.Elems.meanDataAtIndices(x(:,2)) 2*ones(mesh.nElems,1)} ;
%tests(end+1,:) = {sum(x.^2,2) 2*mesh.Elems.meanDataAtIndices(x(:,1)) 2*mesh.Elems.meanDataAtIndices(x(:,2)) 0*mesh.Elems.meanDataAtIndices(x(:,2))} ;
% Apply
errors = NaN(size(tests,1),mesh.nCoord^2) ;
for tt = 1:size(tests,1)
    f = tests{tt,1} ;
    for cc = 1:min(mesh.nCoord^2,size(tests,2)-1)
        d2f_dx2c = G2{cc}*f ;
        errors(tt,cc) = norm(tests{tt,cc+1}-d2f_dx2c) ;
    end
end
errors


