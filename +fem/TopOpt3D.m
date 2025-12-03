%% TOPOLOGY OPTIMISATION FOR COMPLIANCE MINIMISATION
% Highly inspirated from: https://www.topopt.mek.dtu.dk/Apps-and-software/A-99-line-topology-optimization-code-written-in-MATLAB
% Principle: minimize the compliance c(x) = U'*K(x)*U where:
%   - K(x)*U = F
%   - 0 < xmin < x < 1 where x is the material "density" 
%   - int(x)/V0 = f, with V0 the domain volume and f the target volume fraction

%% DEFINE THE BASE GEOMETRY
clc
clear all
import pkg.geometry.levelset.*
    switch 3
        case 1 % 3D cantilever beam
            % Mesh
                L = [2 1 1] ;
                dx = 1/15 ; 
                mesh = pkg.geometry.mesh.GridMesh(L.*[0;1],dx) ;
            % Left face clamped
                fixDOF = mesh.near([0 NaN NaN]).*[1 1 1] ;
            % Loads on right face
                F = mesh.near([max(mesh.Nodes(:,1)) NaN NaN]).*[0 0 1] ;
            % Imposed densities
                xFull = any(fixDOF|F,2) ;
                xZero = false(mesh.nNodes,1) ;
        case 2 % hollow disk
            % Mesh
                Rmin = .3 ; Rmax = 1 ; e = 1/40 ; dR = (Rmax-Rmin)/50 ;
                lvlst = Circle([0 0],Rmax) - Circle([0 0],Rmin) ; 
                mesh = lvlst.mesh(dR) ;
             % Inner radius clamped
                fixDOF = mesh.near(@(p)sum(p.^2,2)-Rmin^2).*[1 1] ;
            % Tangent & Normal load on outer radius
                io = find(mesh.near(@(p)sum(p.^2,2)-Rmax^2)) ;
                xo = mesh.Nodes(io,:) ;
                Fn = sparse(io*[1 1],[1 2].*ones(numel(io),1),xo,mesh.nNodes,2) ;
                Ft = sparse(io*[1 1],[1 2].*ones(numel(io),1),xo*[0 1 ; -1 0],mesh.nNodes,2) ;
            % Harmonic load cases
                if 1
                    nHarm = 4 ; 
                    theta = atan2(mesh.Nodes(:,2),mesh.Nodes(:,1)) ; 
                    modu = cos(nHarm*theta+(0:nHarm-(nHarm>1))*pi/(nHarm+(nHarm==1))).^2 ;
                    F = repmat(Ft,[1 size(modu,2)]).*repelem(modu,1,mesh.nCoord) ;
                    F = [F repmat(Fn,[1 size(modu,2)]).*repelem(modu,1,mesh.nCoord)] ; 
                else
                    Fi = Ft ; 
                    F = [diag(Fi(:,1)) ; diag(Fi(:,2))] ; % one load case for every outer node
                end
                F(:,sum(abs(F),1)==0) = [] ;
            % Cost function
                cFcn = @(c) sum(c(1:end/2))/sum(c(end/2+1:end)) ; 
                dcFcn_dx = @(c,dc_dx) (sum(dc_dx(:,1:end/2),2)*sum(c(end/2+1:end))-sum(c(1:end/2))*sum(dc_dx(:,end/2+1:end),2))/(sum(c(end/2+1:end))^2) ;
            % Imposed densities
                xFull = any(fixDOF,2) | any(reshape(F,mesh.nNodes,[]),2) ...
                            | sum(mesh.Nodes.^2,2)-(Rmin+e)^2<0 ...
                            | sum(mesh.Nodes.^2,2)-(Rmax-e)^2>0 ;
                xZero = false(mesh.nNodes,1) ;
        case 3 % link
            % Mesh
                D = 2 ; R1 = 2/10 ; R2 = 1/10 ; d = 4/10 ; dx = D/30 ; e = 1*dx ; H = .5 ; subdiv = 1 ;
                lvlst = Polygon([0 -R1-d ; 0 R1+d ; D R2+d ; D -R2-d]) ;
                lvlst = lvlst + Circle([0 0],R1+d) + Circle([D 0],R2+d) ;
                lvlst = lvlst - Circle([0 0],R1) - Circle([D 0],R2) ;
                mesh = lvlst.mesh(dx) ;
                mesh.CatmullClark(subdiv) ;
                mesh = mesh.extrude([0 0 H],floor((2^subdiv)*H/dx)) ;
            % First hole is supporting
                fixDOF = mesh.near(@(p)sum(p(:,1:2).^2,2)-R1^2).*[1 1 1] ;
            % Second hole is loaded vertically on its bottom half
                io = find(mesh.near(@(p)(sum((p(:,1:2)-[D 0]).^2,2)-R2^2) + double(p(:,2)>0))) ;
                xo = mesh.Nodes(io,:) ;
                fo = xo.*[0 1 0] ;
                F = sparse(io*[1 1 1],[1 2 3].*ones(numel(io),1),fo,mesh.nNodes,3) ;
            % Imposed densities
                xFull = any(fixDOF|abs(F),2) ...
                            | sum(mesh.Nodes(:,1:2).^2,2)-(R1+e)^2<0 ...
                            | sum((mesh.Nodes(:,1:2)-[D 0]).^2,2)-(R2+e)^2<0 ;
                xZero = false(mesh.nNodes,1) ;
        case 4 % hook
            % Shape params
                H = 17/20 ; 
                R1m = 2/20 ; R1p = 18/100 ;
                R2m = 4/20 ; R2p = 9/20 ; xC2 = 18/100 ;
                e = (R1p-R1m)/3 ; 
                t = 3/10 ;
                dx0 = H/20 ; 
                subdiv = 1 ;
            % Determine remaining parameters
                dd = @(Px,Py,R3m,R3p)[sqrt(Px^2+(H+Py)^2)-(R2m+R3p) ; sqrt((Px-xC2)^2+(H+Py)^2)-(R2p+R3m) ; sqrt(Px^2+Py^2)-(R1p+R3m) ; 2*R1p-(R3p-R3m)] ;
                pp = [H 0 H-R1p H+R1p] ;
                pp = fminsearch(@(x)norm(dd(x(1),x(2),x(3),x(4))),pp,optimset('TolX',1e-9)) ;
                P = pp(1:2) ; R3m = pp(3) ; R3p = pp(4) ;
            % Individual shapes
                clf ; axis equal ; 
                pt = [P ; -(1.1*R1p)/(R1p+R3m)*P.*[1 1 ; 0 1] + [0 0 ; xC2-(1.1*R2p) 0] ; xC2-(1.1*R2p) -H] ;
                poly1 = Polygon([pt ; 0 -H]) ;
                poly2 = Polygon([pt ; xC2 -H]) ;
                C1 = Circle([0 0],R1p) - Circle([0 0],R1m) ;
                C2 = Circle([xC2 -H],R2p) - Circle([0 -H],R2m) ;
                C3 = Circle(P,R3p) - Circle(P,R3m) ;
                %for ll = [poly1 poly2 C1 C2 C3] ; plot(ll) ; end ; return ;
                lvlst = (C3 & poly2) ;
                lvlst = lvlst | (C2 - poly1) ;
                lvlst = (lvlst + Circle([0 0],R1p)) - Circle([0 0],R1m) ;
                lvlst = lvlst + Circle([-(R2p+R2m-xC2)/2 -H],(R2p-R2m-xC2)/2) ;
                %plot(lvlst) ;
                mesh = lvlst.mesh(dx0) ;
                mesh.CatmullClark(subdiv) ;
                mesh.extrude([0 0 t],floor((2^subdiv)*t/dx0))
            % First hole is supporting on its top half
                fixDOF = (mesh.near(@(p)sum(p(:,1:2).^2,2)-R1m^2) & mesh.Nodes(:,2)>=0).*[1 1 1] ;
            % Second hole is loaded vertically on its bottom half
                io = find(mesh.near(@(p)sum((p(:,1:2)-[0 -H]).^2,2)-R2m^2) & mesh.Nodes(:,2)<=-H) ;
                xo = mesh.Nodes(io,:) ;
                fo = xo.*[0 1 0] ;
                F = sparse(io*[1 1 1],[1 2 3].*ones(numel(io),1),fo,mesh.nNodes,3) ;
            % Imposed densities
                xFull = any(fixDOF|F,2) ...
                            | sum(mesh.Nodes(:,1:2).^2,2)-(R1m+e)^2<0 ...
                            | sum((mesh.Nodes(:,1:2)-[0 -H]).^2,2)-(R2m+e)^2<0 & mesh.Nodes(:,2)<=-H ;
                xZero = false(mesh.nNodes,1) ;
        case 5 % cantilever beam
            L = [1 1] ;
            dx = 1/50 ;
            e = 2*dx ;
            % Mesh 
                [xx,yy] = meshgrid(0:dx:L(1),0:dx:L(2)) ;
                mesh = pkg.geometry.mesh.Mesh(cat(3,xx,yy)) ;
            % Left edge clamped
                fixDOF = mesh.near([0 NaN]).*[1 1] ;
            % Loads on right edge
                %F = mesh.near([max(mesh.Nodes(:,1)) NaN]).*cat(3,[0 1],[1 0]) ;
                F = kron(speye(2),diag(sparse(mesh.near([max(mesh.Nodes(:,1)) NaN])))) ;
                F(:,sum(abs(F),1)==0) = [] ;
                F = circshift(F,size(F,2)/2,2) ;
            % Cost function
                cFcn = @(c) sum(c(1:end/2))/sum(c(end/2+1:end)) ; 
                dcFcn_dx = @(c,dc_dx) (sum(dc_dx(:,1:end/2),2)*sum(c(end/2+1:end))-sum(c(1:end/2))*sum(dc_dx(:,end/2+1:end),2))/(sum(c(end/2+1:end))^2) ;
                %cFcn = @(c) exp(sum(c(1:end/2))-sum(c(end/2+1:end))) ; 
                %dcFcn_dx = @(c,dc_dx) (sum(dc_dx(:,1:end/2),2)-sum(dc_dx(:,end/2+1:end),2))*cFcn(c) ;
            % Imposed densities
                xFull = mesh.Nodes(:,1)<e | L(1)-mesh.Nodes(:,1)<e ;
                xZero = false(mesh.nNodes,1) ;
        case 6 % unit cell
            L = [1 1] ;
            dx = 1/80 ;
            %EPS0 = cat(3,0*[1 0 0],0*[0 1 0],[1 1 0],[0 0 1]) ; % macroscopic strain
            EPS0 = cat(3,[0 0 1],[1 1 0]) ; % macroscopic strain
            %EPS0 = cat(3,[1 0 0],[0 1 0]) ; % macroscopic strain
            % Mesh 
                [xx,yy] = meshgrid(0:dx:L(1),0:dx:L(2)) ;
                mesh = pkg.geometry.mesh.Mesh(cat(3,xx,yy)) ;
            % Periodic Boundary Conditions
                Tdof = mesh.perNodeMat ;
                Tdof(:,1) = 0 ;
                Tdof = kron(speye(mesh.nCoord),Tdof) ;
            % Zero forces
                F = zeros(mesh.nNodes,mesh.nCoord,size(EPS0,3)) ;
            % Cost function
                cFcn = @(c) sum(c(:,1:end/2))./sum(c(:,end/2+1:end)) ; 
                dcFcn_dx = @(c,dc_dx) (sum(dc_dx(:,1:end/2),2).*sum(c(:,end/2+1:end),2)-sum(c(:,1:end/2),2)*sum(dc_dx(:,end/2+1:end),2))./(sum(c(:,end/2+1:end))^2) ;
%                 n = 2 ; cFcn = @(c) c(:,1).^n + c(:,2).^-n ; dcFcn_dx = @(c,dc_dx) n*(dc_dx(:,1).*c(:,1).^(n-1)-dc_dx(:,2).*c(:,2).^(-n-1)) ;
            % Imposed densities
                xFull = false(mesh.nNodes,1) ;
                xZero = false(mesh.nNodes,1) ;
        otherwise % template
            mesh = pkg.geometry.mesh.Mesh() ; % the mesh
            dofSZ = [mesh.nNodes mesh.nCoord] ;
            fixDOF = false(dofSZ) ; % fixed DOFs
            F = zeros(dofSZ) ; % applied load
            xFull = false(dofSZ) ; % density imposed to 1
            xZero = false(dofSZ) ; % density imposed to 0
    end
% Other information
    if ~exist('Tdof','var') ; Tdof = diag(sparse(~fixDOF(:))) ; end
    Tdof(:,sum(Tdof,1)==0) = [] ;
    nLoadCases = numel(F)/mesh.nNodes/mesh.nCoord ;
% Display mesh
    clf ; axis equal ; pl = plot(mesh,'VisibleFaces','all') ;
    if exist('fixDOF','var') ; pl.Selected.Nodes = find(any(fixDOF,2)) ; end
    pl.CData = [1 1 1] - [1 0 1].*xFull ;
% Display loads
    Fq = reshape(F,[],nLoadCases) ; 
    [iq,jq] = find(repmat(speye(mesh.nNodes),[1 mesh.nCoord])*Fq) ;
    iiq = sub2ind(size(Fq),iq+(0:mesh.nCoord-1)*mesh.nNodes,repmat(jq,[1 mesh.nCoord])) ;
    Fq = Fq(iiq) ; Xq = mesh.Nodes(iq,:) ;
    pl.Highlighted.Nodes = iq ;
    q = quiver(Xq(:,1),Xq(:,2),Fq(:,1),Fq(:,2),'b','linewidth',1.5) ;

%% TOPOLOGY OPTIMIZATION
    nu = 0.3 ; % Poisson coefficient
    penal = 3 ; % density penalization power d = x^p
    volfrac = .25 ; % volume fraction 
    xmin = 1e-3 ; % minimum material density
    initialDistribution = @(N)ones(N)+0.1*rand(N) ; % initial material distribution
    maxChange = 0.02 ; % maximum change in material density by iteration
    eta = .1 ; % damping factor
    smoothRatio = .15 ; % weight of neightbors
    smoothPower = 2 ; % Recursive smoothing
    minChange = 0.001 ; % minimum change in density at convergence
    maxIt = 500 ; % maximum number of iterations
    plotFreq = 1 ; % plot update frequency
% Integration matrices
    [ee,w,ie] = mesh.integration ;
    V0 = sum(w) ; % mesh volume
    nGP = numel(w) ;
    N = mesh.interpMat(ee,ie) ;
    G = mesh.diffMat(ee,ie) ;
    O = sparse(nGP,mesh.nNodes) ;
% Strain matrix [E11;E22;2E12] = B*u(:)
    B = [G{1} O O ; O G{2} O ; O O G{3} ; G{2} G{1} O ; G{3} O G{1} ; O G{3} G{2}] ; 
% Plane stiffness
    %C = sparse([1-nu nu 0 ; nu 1-nu 0 ; 0 0 (1-2*nu)/2]) ; % plane strain
    % C = sparse([1 nu 0 ; nu 1 0 ; 0 0 (1-nu)/2]) ; % plane stress
    C = pkg.fem.bloch.stiffness(1,.5/(1+nu)) ;
% Boundary conditions
    Bk = B*Tdof ;
    Fk = Tdof'*reshape(F,[],nLoadCases) ;
% Density smoothing filter
    S = sparse(mesh.Edges.NodeIdx(:,[1 2]),mesh.Edges.NodeIdx(:,[2 1]),1,mesh.nNodes,mesh.nNodes) ;
    S = speye(mesh.nNodes) + smoothRatio*S./sum(S,2)  ;
    S = S./sum(S,2) ;
    S = S^smoothPower ;
    Nx = N*S ;
% Volume derivative
    dv_dx = sum(Nx'*diag(sparse(w)),2) ;
% Initialization
    x = volfrac.*initialDistribution([size(Nx,2),1]) ;
    x(xFull) = 1 ;
    x(xZero) = xmin ;
    x = S*x ;
    U = zeros(mesh.nNodes,mesh.nCoord,nLoadCases) ;
    clf ; axis equal tight ; pl =  plot(mesh,'VisibleEdges','none','VisibleFaces','all') ; view(35,30)
% UPDATE LOOP
    profile on
    it = 0 ; change = Inf ; outFlag = '' ; lastPlotTime = tic ;
    while isempty(outFlag)
    % FE stiffness matrix
        xe = Nx*x ;
        W = kron(C,diag(sparse(w.*xe.^penal))) ;
        K = Bk'*W*Bk ;
        K = K+K' ;
    % Force vector
        f = Fk ;
        if exist('EPS0','var')
%             f = f - Bk'*W*kron(speye(3),N)*reshape(EPS0.*ones(mesh.nNodes,1),[],nLoadCases) ;
            f = f - Bk'*W*reshape(EPS0.*ones(nGP,1),[],nLoadCases) ;
        end
    % Displacement
        Uk = K\f ;
        U = reshape(full(Tdof*Uk),mesh.nNodes,mesh.nCoord,nLoadCases) ;
    % Strains & Stresses
        EPS = reshape(B*reshape(U,[],nLoadCases),nGP,size(C,1),nLoadCases) ;
        if exist('EPS0','var') ; EPS = EPS+EPS0 ; end
        SIG = (xe.^penal).*reshape(sum(reshape(full(C),[1 size(C)]).*reshape(EPS,nGP,1,size(C,1),nLoadCases),3),nGP,size(C,1),nLoadCases) ;
    % Compliance & derivative 
        ce = sum(w.*SIG.*EPS.*[1 1 1 2 2 2],2) ;
        c = reshape(sum(ce,1),[1 nLoadCases]) ;% sum((K*Uk).*Uk,1)' ;%
        %dce_dx = reshape(kron(speye(nLoadCases),Nx')*diag(sparse(reshape((penal./xe).*ce,[],1))),[],nLoadCases) ;
        dc_dx = Nx'*reshape((penal./xe).*ce,[nGP nLoadCases]) ;% Nx'*(penal*xe.^(penal-1).*w.*reshape(sum(EPS.^2,2),[],nLoadCases)) ;%
    % Cost function derivative
        if exist('cFcn','var')
            fcn = cFcn(c) ;
            dfcn_dx = dcFcn_dx(c,dc_dx) ;
%             fcne = cFcn(ce(:,:)) ;
%             dfcne_dx = reshape(dcFcn_dx(repelem(ce(:,:),mesh.nNodes,1),dce_dx),[],nGP) ;
        else % sum of compliances corresponding to load cases
            fcn = sum(c) ;
            dfcn_dx = sum(dc_dx,2) ;
%             fcne = sum(ce(:,:),2) ;
%             dfcne_dx = reshape(sum(dce_dx,2),[],nGP) ;
        end
    % Optimal Criteria update
        l1 = 0; l2 = 1e9 ; xnew = x ; 
        d_dx = dfcn_dx./dv_dx ;
        xBen = x.*d_dx.*abs(d_dx).^(eta-1) ;
        iii = 0 ;
        while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1);
            xnew = xBen/lmid^eta ;
            xnew = min(x+maxChange,xnew) ;
            xnew = min(1,xnew) ;
            xnew = max(x-maxChange,xnew) ;
            xnew = max(xmin,xnew);
            xnew(xFull) = 1 ;
            xnew(xZero) = xmin ;
            V = sum((Nx*xnew).*w) ;
            if V > volfrac*V0 ; l1 = lmid ; else l2 = lmid ; end
            iii=iii+1;
        end
        change = max(abs(xnew(:)-x(:)));
        x = S*xnew;%
    % Convergence ?
        it = it+1 ;
        if change<minChange ; outFlag = 'minChange' ; end
        if it>=maxIt ; outFlag = 'maxIt' ; end
    % Display
        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',it,full(fcn),full(V),full(change));
        if ~isempty(outFlag) || toc(lastPlotTime)>1/plotFreq
            %pl.Deformation = 0.1*U/max(abs(U(:))).*norm(range(mesh.boundingBox,1)) ;
            VM = sqrt(sum(SIG.^2,2)) ;
            pl.CData = sum(VM(:,:,end),3) ;
            pl.Faces.FaceAlpha = 'interp' ;
            pl.Faces.FaceVertexAlphaData = full(x) ;
            %pl.CData = x ;
            drawnow ;
            lastPlotTime = tic ;
        end
    end
    disp(outFlag) ;
    profile off
    return
    
%% PLOT THE SOLUTION
clf ; axis equal ; plot(mesh,'Deformation',U(:,:,end)) ;
    
%% REPLICATE THE MESH
N = [1 1]*5 ; 
dx = range(mesh.Nodes,1) ;
[xx,yy] = meshgrid((0:N(1)-1)*dx(1),(0:N(2)-1)*dx(2)) ;
mesh.move([xx(:) yy(:)]) ;
[mesh,~,~,meanMat] = mesh.cullDuplicates() ;
x = meanMat*repmat(x,prod(N),1) ;

clf ; axis equal ; plot(mesh,'CData',x,'VisibleEdges','none') ;

%% CUT THE SHAPE
    lvl = 30/100 ;
% Extract inner boundaries
   shape = mesh.simplex().cut(lvl-x).IN ;
   shape.clean ;
% Display
    clf ; axis equal ; pl = plot(shape,'EdgeColor','none') ;
    light
    view(35,30) ;
    
%% CUT AND SMOOTH THE SHAPE
    lvl = 30/100 ;
    smoothIt = 100 ;
% Extract inner boundaries
    shape = mesh.simplex().cut(lvl-x).IN ;
    shape.Elems = shape.Faces.subpart(shape.outerFaces) ;
% Laplacian smoothing
    fixNodes = mesh.Nodes(false(mesh.nNodes,1) | xFull,:) ;
    [pp,d2] = shape.closestNode(fixNodes) ;
    pp = pp(d2<shape.defaultTolerance) ;
    smooth = ~full(logical(sparse(pp,pp*0+1,1,shape.nNodes,1))) ;% & ~xZero ;
    e2n = shape.elem2node ;
    me = e2n*diag(sparse(1./sum(e2n,1))) ;
    mn = diag(sparse(1./sum(e2n,2)))*e2n ;
    mn = mn(smooth,:) ;
    for ii = 1:smoothIt
        xe = me'*shape.Nodes ; % elem 
        shape.Nodes(smooth,:) = (shape.Nodes(smooth,:) + mn*xe)/2 ;
        ii
    end
% % Add original mesh's shape
%     out = pkg.geometry.mesh.Mesh('Nodes',mesh.Nodes,'Elems',mesh.Edges.subpart(mesh.boundaryEdges)) ;
%     out.Nodes(x<lvl,:) = NaN ;
%     out.clean()
%     shape = shape.merge(out).clean ;
% % Regenerate the 2D mesh from boundaries (delaunay triang.)
%     shape = pkg.geometry.mesh.Mesh('Nodes',shape.Nodes,'Elems',delaunay(shape.Nodes)) ;
%     valid = mesh.interpMat(shape.centroid,[],mesh.Elems,false)*x>lvl ;
%     shape.Elems = shape.Elems.subpart(valid) ;
% Display
    clf ; axis equal ; pl = plot(shape,'EdgeColor','none') ;
    light
    view(35,30) ;
    
%% CUT AND SMOOTH THE SHAPE
    lvl = 40/100 ;
    smoothIt = 50 ;
% Extract inner boundaries
    shape = mesh.simplex().cut(lvl-x).ON ;
    shape.clean ;
% Laplacian smoothing
    smooth = true(shape.nNodes,1) & ~shape.outerNodes ;
    e2n = shape.elem2node ;
    me = e2n./sum(e2n,1) ;
    mn = e2n./sum(e2n,2) ;
    for ii = 1:smoothIt
        xe = me'*shape.Nodes ; % elem 
        shape.Nodes(smooth,:) = (shape.Nodes(smooth,:) + mn(smooth,:)*xe)/2 ;
    end
% % Add original mesh's shape
%     out = pkg.geometry.mesh.Mesh('Nodes',mesh.Nodes,'Elems',mesh.Edges.subpart(mesh.boundaryEdges)) ;
%     out.Nodes(x<lvl,:) = NaN ;
%     out.clean()
%     shape = shape.merge(out).clean ;
% % Regenerate the 2D mesh from boundaries (delaunay triang.)
%     shape = pkg.geometry.mesh.Mesh('Nodes',shape.Nodes,'Elems',delaunay(shape.Nodes)) ;
%     valid = mesh.interpMat(shape.centroid,[],mesh.Elems,false)*x>lvl ;
%     shape.Elems = shape.Elems.subpart(valid) ;
% Display
    clf ; axis equal ; pl = plot(shape,'EdgeColor','none') ;
    light
    view(35,30) ;
    
    
%% BUILD THE PART
    thickness = 5/50 ;
    plateThickness = 0/50 ;
    scale = 50 ;
    part = shape.extrude([0 0 thickness]) ;
    if plateThickness ; part = part.move([0 0 plateThickness]).merge(mesh.extrude([0 0 plateThickness])) ; end
    part.Nodes = part.Nodes*scale ;
    clf ; axis equal ; plot(part,'VisibleEdges','none')%,'FaceColor',[1 1 1]*.7) ;
    %light() ; lighting gouraud
    
%% EXPORT AS STL
    part.stlwrite() ;
    


