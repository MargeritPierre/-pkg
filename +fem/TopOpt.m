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
    switch 4
        case 1 % cantilever beam
        % Mesh
            dx = 1/50 ; 
            [xx,yy] = meshgrid(0:dx:1,0:dx:1) ;
            mesh = pkg.geometry.mesh.Mesh(cat(3,xx,yy)) ;
        % Left edge clamped
            fixDOF = mesh.near([0 NaN]).*[1 1] ;
        % Loads on right edge
            F = mesh.near([max(mesh.Nodes(:,1)) NaN]).*[0 -1] ;
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
        % Tangent load on outer radius
            io = find(mesh.near(@(p)sum(p.^2,2)-Rmax^2)) ;
            xo = mesh.Nodes(io,:) ;
            fo = xo*(0*eye(2) + [0 1 ; -1 0]) ;
            F = sparse(io*[1 1],[1 2].*ones(numel(io),1),fo,mesh.nNodes,2) ;
        % Imposed densities
            xFull = any(fixDOF|F,2) ...
                        | sum(mesh.Nodes.^2,2)-(Rmin+e)^2<0 ...
                        | sum(mesh.Nodes.^2,2)-(Rmax-e)^2>0 ;
            xZero = false(mesh.nNodes,1) ;
        case 3 % link
        % Mesh
            D = 2 ; R1 = 2/10 ; R2 = 1/10 ; d = 1/20 ; dx = D/200 ; e = 3*dx ; 1/40 ;
            lvlst = Polygon([0 -R1-d ; 0 R1+d ; D R2+d ; D -R2-d]) ;
            lvlst = lvlst + Circle([0 0],R1+d) + Circle([D 0],R2+d) ;
            lvlst = lvlst - Circle([0 0],R1) - Circle([D 0],R2) ;
            mesh = lvlst.mesh(dx) ;
        % First hole is supporting
            fixDOF = mesh.near(@(p)sum(p.^2,2)-R1^2).*[1 1] ;
        % Second hole is loaded vertically on its bottom half
            io = find(mesh.near(@(p)(sum((p-[D 0]).^2,2)-R2^2) + double(p(:,2)>0))) ;
            xo = mesh.Nodes(io,:) ;
            fo = xo.*[0 1] ;
            F = sparse(io*[1 1],[1 2].*ones(numel(io),1),fo,mesh.nNodes,2) ;
        % Imposed densities
            xFull = any(fixDOF|F,2) ...
                        | sum(mesh.Nodes.^2,2)-(R1+e)^2<0 ...
                        | sum((mesh.Nodes-[D 0]).^2,2)-(R2+e)^2<0 ;
            xZero = false(mesh.nNodes,1) ;
        case 4 % hook
        % Shape params
            H = 17/20 ; 
            R1m = 2/20 ; R1p = 18/100 ;
            R2m = 4/20 ; R2p = 9/20 ; xC2 = 18/100 ;
            e = (R1p-R1m)/3 ; 
            dx = H/100 ;
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
            mesh = lvlst.mesh(dx) ;
            plot(mesh) ;
        % First hole is supporting on its top half
            fixDOF = (mesh.near(@(p)sum(p.^2,2)-R1m^2) & mesh.Nodes(:,2)>=0).*[1 1] ;
        % Second hole is loaded vertically on its bottom half
            io = find(mesh.near(@(p)sum((p-[0 -H]).^2,2)-R2m^2) & mesh.Nodes(:,2)<=-H) ;
            xo = mesh.Nodes(io,:) ;
            fo = xo.*[0 1] ;
            F = sparse(io*[1 1],[1 2].*ones(numel(io),1),fo,mesh.nNodes,2) ;
        % Imposed densities
            xFull = any(fixDOF|F,2) ...
                        | sum(mesh.Nodes.^2,2)-(R1m+e)^2<0 ...
                        | sum((mesh.Nodes-[0 -H]).^2,2)-(R2m+e)^2<0 & mesh.Nodes(:,2)<=-H ;
            xZero = false(mesh.nNodes,1) ;
    end
    clf ; axis equal ; pl = plot(mesh) ;
    pl.Selected.Nodes = find(any(fixDOF,2)) ;
    pl.Highlighted.Nodes = find(any(F,2)) ;
    pl.CData = [1 1 1] - [1 0 1].*xFull ;
    quiver(mesh.Nodes(:,1),mesh.Nodes(:,2),F(:,1),F(:,2),'b','linewidth',1.5)%,'AutoScale','off') ;

%% TOPOLOGY OPTIMIZATION
    nu = 0.3 ; % Poisson coefficient
    penal = 3 ; % density penalization power d = x^p
    volfrac = 0.25 ; % volume fraction 
    xmin = 1e-6 ; % minimum material density
    maxChange = 0.2 ; % maximum change in material density by iteration
    eta = .5 ; % damping factor
    minChange = 0.01 ; % minimum change in density at convergence
    maxIt = 200 ; % maximum number of iterations
    plotFreq = 3 ; % plot update frequency
% Integration matrices
    [ee,w,ie] = mesh.integration ;
    V0 = sum(w) ; % mesh volume
    nGP = numel(w) ;
    N = mesh.interpMat(ee,ie) ;
    G = mesh.diffMat(ee,ie) ;
    O = sparse(nGP,mesh.nNodes) ;
% Strain matrix [E11;E22;2E12] = B*u(:)
    B = [G{1} O ; O G{2} ; G{2} G{1}] ; 
% Plane stress stiffness
    C = sparse([1-nu nu 0 ; nu 1-nu 0 ; 0 0 1-2*nu]) ;
% Boundary conditions
    keepDOF = ~fixDOF ;
    Bk = B(:,keepDOF) ;
    Fk = F(keepDOF) ;
% Density smoothing filter
    S = sparse(mesh.Edges.NodeIdx(:,[1 2]),mesh.Edges.NodeIdx(:,[2 1]),1,mesh.nNodes,mesh.nNodes) ;
    S = speye(mesh.nNodes) + 0.9*S./sum(S,2)  ;
    S = S./sum(S,2) ;
    S = S^1 ;
    N = N*S' ;
% Volume derivative
    dv_dx = sum(N'*diag(sparse(w)),2) ; ones(size(N,2),1) ;
% Initialization
    x = volfrac.*ones(size(N,2),1) ;
    x(xFull) = 1 ;
    x(xZero) = xmin ;
    x = S*x ;
    U = zeros(mesh.nNodes,2) ;
    clf ; axis equal ; pl =  plot(mesh,'VisibleEdges','none') ;
% UPDATE LOOP
    it = 0 ; change = Inf ; outFlag = '' ; lastPlotTime = tic ;
    while isempty(outFlag)
    % FE stiffness matrix
        xe = N*x ;
        W = kron(C,diag(sparse(w.*xe.^penal))) ;
        K = Bk'*W*Bk ;
        K = K+K' ;
    % Displacement
        Uk = K\Fk ;
        U(keepDOF) = Uk ;
    % Compliance derivative 
        EPS = reshape(B*U(:),nGP,3) ;
        c = Uk'*K*Uk ;
        dc_dx = N'*(penal*xe.^(penal-1).*w.*sum(EPS.^2,2)) ;
    % Optimal Criteria update
        l1 = 0; l2 = 1e9 ; xnew = x ; 
        xBen = x.*(dc_dx./dv_dx).^eta ;
        while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1);
            xnew = xBen/lmid^eta ;
            xnew = min(x+maxChange,xnew) ;
            xnew = min(1,xnew) ;
            xnew = max(x-maxChange,xnew) ;
            xnew = max(xmin,xnew);
            xnew(xFull) = 1 ;
            xnew(xZero) = xmin ;
            V = sum((N*xnew).*w) ;
            if V > volfrac*V0 ; l1 = lmid ; else l2 = lmid ; end
        end
        change = max(abs(xnew(:)-x(:)));
        x = S*xnew;
    % Convergence ?
        it = it+1 ;
        if change<minChange ; outFlag = 'minChange' ; end
        if it>=maxIt ; outFlag = 'maxIt' ; end
    % Display
        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',it,full(c),V,change);
        if ~isempty(outFlag) || toc(lastPlotTime)>1/plotFreq
            %pl.Deformation = 0.1*U/max(abs(U(:))).*norm(range(mesh.boundingBox,1)) ;
            pl.CData = (sqrt(sum((xe.*EPS*C).^2,2))) ;
            pl.Faces.FaceAlpha = 'interp' ;
            pl.Faces.FaceVertexAlphaData = x ;
            %pl.CData = x ;
            drawnow ;
            lastPlotTime = tic ;
        end
    end
    disp(outFlag) ;
    
%% CUT THE SHAPE
    lvl = 20/100 ;
% Extract inner boundaries
    shape = mesh.cut(lvl-x).IN ;
    shape.clean ;
% Display
    clf ; axis equal ; pl = plot(shape,'VisibleEdges','none') ;
    
    
%% CUT AND SMOOTH THE SHAPE
    lvl = 20/100 ;
    smoothIt = 50 ;
% Extract inner boundaries
    shape = mesh.cut(lvl-x).ON ;
    shape.clean ;
% Laplacian smoothing
    e2n = shape.elem2node ;
    me = e2n./sum(e2n,1) ;
    mn = e2n./sum(e2n,2) ;
    for ii = 1:smoothIt
        xe = me'*shape.Nodes ; % elem 
        shape.Nodes = (shape.Nodes + mn*xe)/2 ;
    end
% Add original mesh's shape
    out = pkg.geometry.mesh.Mesh('Nodes',mesh.Nodes,'Elems',mesh.Edges.subpart(mesh.boundaryEdges)) ;
    shape = shape.merge(out).clean ;
% Regenerate the 2D mesh from boundaries (delaunay triang.)
    shape = pkg.geometry.mesh.Mesh('Nodes',shape.Nodes,'Elems',delaunay(shape.Nodes)) ;
    valid = mesh.interpMat(shape.centroid,[],mesh.Elems,false)*x>lvl ;
    shape.Elems = shape.Elems.subpart(valid) ;
% Display
    clf ; axis equal ; pl = plot(shape,'VisibleEdges','none') ;
    
    
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
    


