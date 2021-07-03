classdef Mesh < pkg.geometry.levelset.LevelSet
% Level set built from a mesh

properties (SetAccess=immutable)
    LevelSetMesh pkg.geometry.mesh.Mesh
end

properties
    LevelSetValues (:,1)
    AutoBuildEdges (1,1) logical = true ;
end

properties (AbortSet)
    Thickness (1,1) = 0 ;
    Theta (1,1) = .5 ;
    AutoReInit (1,1) logical = true ;
end

methods
    function this = Mesh(mesh,lvlstval,varargin)
    % Class construtor
        this = this@pkg.geometry.levelset.LevelSet() ;
        if nargin<2 ; error('Not enougth arguments provided.') ; end
        this.LevelSetMesh = mesh ;
    % Set common levelset properties
        this.BoundingBox = this.LevelSetMesh.boundingBox ;
        this.Kinks = [] ;
        this.EdgeFcns = {} ;
    % Set custom properties
        for arg = 1:2:numel(varargin)-1
            this.(varargin{arg}) = varargin{arg+1} ;
        end
    % Set the level set values 
    %   - set distance function
    %   - may build the boundaries edges
        this.LevelSetValues = lvlstval ;
    end
    
    function d = distanceFun(this,p)
    % Signed distance function interpolation
        N = this.LevelSetMesh.interpMat(p,[],this.LevelSetMesh.Elems,true) ;
        d = N*this.LevelSetValues ;
    end
    
    function [this,edgP] = buildEdges(this)
    % Build the levelset edges using Mesh slice
    % edgP contains all edge points sorted as curves, with NaNs splitting
    % the different curves
        lvl = this.LevelSetValues ;
        if this.isCharacteristic ; lvl = 0.5-lvl ; end
        LVLST = this.LevelSetMesh.cut(lvl).IN ;
    % Separate the mesh curves
        [~,crvs] = LVLST.boundaryCurves ;
    % Convert to edge functions
        this.EdgeFcns = {} ;
        points = {} ;
        for ee = 1:numel(crvs)
            points{ee} = LVLST.Nodes(crvs{ee},:) ;
            % clean the point list (because of the bad mesh.cut function)
                [~,ind] = uniquetol(points{ee},this.LevelSetMesh.defaultTolerance,'byrows',true,'datascale',1) ;
                points{ee} = points{ee}(sort(ind),:) ; 
                points{ee} = points{ee}([1:end,1],:) ; 
            this.EdgeFcns{end+1} = pkg.geometry.levelset.polylineEdgeFunctions(points{ee},'full') ;
        end
    % Return edge points if needed
        if nargout<2 ; return ; end
        points = [points(:)' ; repmat({NaN*(1:LVLST.nCoord)},1,numel(points))] ;
        edgP = cat(1,points{:}) ;
    end
    
    function h = plot(this,varargin)
    % Override the plot function to plot teh level set as mesh
        isAsMesh = cellfun(@(arg)ischar(arg) && strcmpi(arg,'asmesh'),varargin) ;
        if any(isAsMesh)
            varargin(isAsMesh) = [] ;
            if this.LevelSetMesh.nCoord<3
                varargin(end+(1:2)) = {'Deformation',[zeros(size(this.LevelSetMesh.Nodes)) this.LevelSetValues]} ;
            else
                varargin(end+(1:2)) = {'CData',this.LevelSetValues} ;
            end
            h = this.LevelSetMesh.plot(varargin{:}) ;
        else
            h = plot@pkg.geometry.levelset.LevelSet(this,varargin{:}) ;
        end
    end
end

%% LEVEL SET ADVECTION
methods
    function [this,keepDOF] = advect(this,U)
    % Advection of the levelset with the nodal displacement U
    % Solve dphi = - U.grad(phi) using a variationnal formulation
    % with s virtual field: 
    %   int(s*dphi*dx) =
    %       int(s*div(U)*phi*dx) 
    %       - int(s*(U.n)*phi*dS) 
    %       + int(grad(s).u*phi.dx)
    % with dphi = phi_{n+1} - phi{n} and phi = theta*phi_{n+1} + (1-theta)*phi{n}
    % keepDOF [nNodes 1] logical, imposed BCs
    % Domain integrals
        [ee,we,ie] = this.LevelSetMesh.integration ;
        W = diag(sparse(we)) ;
        N = this.LevelSetMesh.interpMat(ee,ie) ;
        G = this.LevelSetMesh.diffMat(ee,ie) ;
    % int(s*dphi*dx)
        Ksp = N'*W*N ;
    % int(s*div(U)*phi*dx)
        divU = cat(2,G{:})*U(:) ;
        Kdiv = N'*diag(sparse(we.*divU))*N ;
    % int(grad(s).u*phi.dx)
        u = N*U ;
        Kgsup = cat(1,G{:})'*diag(sparse(u(:)))*repmat(W,numel(G),1)*N ;
    % Boundary integral: int(s*(U.n)*phi*dS)
        [eb,wb,ib,normal] = this.LevelSetMesh.boundaryIntegration ;
        Nb = this.LevelSetMesh.interpMat(eb,ib) ;
        un = sum((Nb*U).*normal,2) ;
        Ksunp = Nb'*diag(sparse(wb.*un))*Nb ;
    % Dirichlet BC on the inflow boundaries (U.n<0 or U==0)
        isInflow = logical(Nb'*(un<0)) ;
        fixDOF = isInflow ; % | all(U==0,2) ;
        keepDOF = ~fixDOF ;
    % Advection
        dR = Kdiv + Kgsup - Ksunp ; 
        K = Ksp - this.Theta*dR ;
        dphi = K(keepDOF,keepDOF)\(dR(keepDOF,:)*this.LevelSetValues) ;
        this.LevelSetValues(keepDOF) = this.LevelSetValues(keepDOF) + dphi ;
    % Reinitialization ?
        if this.AutoReInit ; this = this.reinit(keepDOF) ; end
    end
end

%% LEVEL SET REINITIALISATION
methods
    function [d,beta] = regularizationParameter(this)
    % Regularization used in the computation of optimal reinit. thickness & time step
        d = .1 ; beta = .5 ;
    end
    
    function dx = deltaX(this)
    % Characteristic length corresponding to spatial discretization
    % Use the median edge length
        dx = median(this.LevelSetMesh.elemSize(this.LevelSetMesh.Edges)) ;
    end
    
    function th = defaultThickness(this)
    % Default interface thickness using the mesh edge length
        [d,beta] = this.regularizationParameter ;
        dx = this.deltaX ;
        th = beta*dx^(1-d) ;
    end
    
    function dt = defaultReinitTimeStep(this)
    % Default virtual time step used in the reinitialization procedure
        [d,beta] = this.regularizationParameter ;
        dx = this.deltaX ;
        dt = beta*dx^(1+d) ;
    end
    
    function this = reinit(this,keepDOF,tol,maxIt,dt)
    % Reinitialization of the levelset up to a tolerance
    % keepDOF is used to impose Dirichlet BCs
    % Characteristic functions only
    % Solve dphi = e*div(grad(phi)) - div(phi*(1-phi).n0)
    % with s virtual field:
    %   int(s*dphi*dx) =
    %       int(grad(s).n0*phi*(1-phi).dx) 
    %       - int(e*grad(s).grad(phi)*dx) 
    %       + int(e*s*n.grad(phi)*dS) 
    %       - int(s*n0.n*phi*(1-phi)*dS) 
    % with dphi = phi_{n+1} - phi{n} and phi = theta*phi_{n+1} + (1-theta)*phi{n}
    % Convert to characteristic function
        isSignedDistance = ~isCharacteristic(this) ;
        if isSignedDistance ; this.Thickness = this.defaultThickness ; end
        epsN = 1e-3/this.Thickness ; % normal regulariation when the gradient vanishes
    % Diriclet BC by default
        if nargin<2
            %keepDOF = ~this.LevelSetMesh.boundaryNodes ; 
            keepDOF = true(this.LevelSetMesh.nNodes,1) ; 
        end
    % Iteration infos
        if nargin<3 ; tol = 1e-2 ; end
        if nargin<4 ; maxIt = 10 ; end
        if nargin<5 ; dt = this.defaultReinitTimeStep ; end
        nCoord = this.LevelSetMesh.nCoord ;
    % Initialization
        phi = this.LevelSetValues ;
        phi = min(max(phi,0),1) ;
    % Domain integrals
        [ee,we,ie] = this.LevelSetMesh.integration ;
        W = diag(sparse(we)) ;
        N = this.LevelSetMesh.interpMat(ee,ie) ;
        G = this.LevelSetMesh.diffMat(ee,ie) ;
        grad = cat(1,G{:}) ;
    % Levelset normal n0 = grad(phi)/|grad(phi)| (on the volume)
        gradPhi = reshape(grad*phi,numel(we),nCoord) ; % [nGP nCoord]
        n0 = gradPhi./(sqrt(sum(gradPhi.^2,2))+epsN) ;
    % int(s*dphi*dx)
        Ksp = N'*W*N ;
        Ksp = diag(sum(Ksp,2)) ;
    % int(grad(s).n0*phi*(1-phi).dx) 
        dRgsn0pp = (diag(sparse(n0(:)))*grad)'*repmat(W,[nCoord 1])*N ;
    % int(e*grad(s).grad(phi)*dx)
        Kgsgp = this.Thickness*grad'*kron(speye(nCoord),W)*grad ;
    % Boundary integrals
        [eb,wb,ib,normal] = this.LevelSetMesh.boundaryIntegration ;
        Wb = diag(sparse(wb)) ;
        Nb = this.LevelSetMesh.interpMat(eb,ib) ;
        Gb = this.LevelSetMesh.diffMat(eb,ib) ;
        gradB = cat(1,Gb{:}) ;
    % Levelset normal n0 = grad(phi)/|grad(phi)| (on the boundary)
        gradPhiB = reshape(gradB*phi,numel(wb),nCoord) ; % [nGP nCoord]
        n0B = gradPhiB./sqrt(sum(gradPhiB.^2,2)+epsN) ;
    % int(e*s*n.grad(phi)*dS) 
        Ksngp = this.Thickness*Nb'*repmat(Wb,[1 nCoord])*diag(sparse(normal(:)))*gradB ;
    % int(s*n0.n*phi*(1-phi)*dS) 
        dRsn0npp = Nb'*diag(sparse(sum(n0B.*normal,2).*wb))*Nb ;
    % Theta-Scheme: M*dphi = r(phi) + theta*dR_dphi(phi)*dphi
        M = (1/dt)*Ksp ; 
        dRpp = dRgsn0pp - dRsn0npp ; % non-linear part
        dRp = - Kgsgp + 0*Ksngp ; % linear part
        dphi = Inf ; it = 0 ;
        while it<maxIt && norm(dphi)>tol
            r = dRpp*(phi.*(1-phi)) + dRp*phi ;
            dR_dphi = dRpp*diag(sparse(1-2*phi)) + dRp ;
            K = M - this.Theta*dR_dphi ;
            dphi = K(keepDOF,keepDOF)\r(keepDOF) ;
            phi(keepDOF) = phi(keepDOF) + dphi ;
            it = it+1 ;
            %cla ; axis equal ; pl = plot(this.LevelSetMesh,'Deformation',[0 0 1].*phi) ; pl.Selected.Nodes = ~keepDOF ; drawnow ;
        end
        %[it norm(dphi) norm(r)]
    % Set the levelset
        this.LevelSetValues = phi ;
    % Back to a signed distance function
        if isSignedDistance ; this.Thickness = 0 ; end
    end
end

%% PROPERTIES SET METHODS
methods
    function this = set.LevelSetValues(this,val)
    % Change the levelSet values
        if numel(val)~=this.LevelSetMesh.nNodes ; error('The level set values are defined on the mesh nodes.') ; end
        this.LevelSetValues = val(:) ;
        this.Function = @this.distanceFun ;
        if this.AutoBuildEdges ; this = this.buildEdges ; end
    end
    
    function this = set.AutoBuildEdges(this,val)
    % Auto-build levelset edge functions on level set update ?
        this.AutoBuildEdges = logical(val) ;
        if this.AutoBuildEdges ; this = this.buildEdges ; end
    end
    
    function this = set.Thickness(this,t)
    % Change the levelSet values
        lvl = this.LevelSetValues ;
        if this.isCharacteristic % back to signed distance function
            lvl = this.Thickness*log((1-lvl)./lvl) ;
        end
        if t~=0 % to characteristic function
            lvl = 1./(1+exp(lvl/t)) ; 
        end 
        this.Thickness = t ;
        this.LevelSetValues = lvl ;
    end
    
    function val = isCharacteristic(this)
    % Is it a charcteristic level-set ?
        val = this.Thickness~=0 ;
    end
end
    
end







%% UNIT TESTS
function tests
%% ROTATING SHAPES
    clc ; %clearvars
    % Parameters
        shape = 'diskslit' ; % see below
        margin = 1/3 ; de = 1/100 ;
        elemType = 'quad' ; % 'tri' or 'quad'
        angle = 2*pi ; nIt = 200 ;
    % Geometry levelset
        import pkg.geometry.levelset.*
        switch lower(shape)
            case 'square'
                C = 1 ; geo = Rectangle([0 0],[C C]) ;
            case 'disk'
                R = 1 ; geo = Circle([0 0],R) ;
            case 'diskslit'
                R = 1 ; l = R/3 ; 
                geo = Circle([0 0],R) - Rectangle([-l/2 -R-margin; l/2 0]) ;
            case 'dumbell'
                D = 1 ; R = D/4 ; H = R/1 ;
                geo = Circle([0 0],R) + Circle([D 0],R) + Rectangle([0 -H/2 ; D H/2]) ;
        end
    % Enhance geometry bounding box
        bbox = geo.BoundingBox ;
        bbox = mean(bbox,1) + (0.5+margin)*[-1;1]*max(range(bbox,1)) ;
        dx = de*min(range(bbox,1)) ;
    % Build the mesh
        switch elemType
            case 'tri'
                mesh = pkg.geometry.levelset.Rectangle(bbox).mesh(dx) ;
            case 'quad'
                [XX,YY] = meshgrid(bbox(1,1):dx:bbox(2,1),bbox(1,2):dx:bbox(2,2)) ;
                mesh = pkg.geometry.mesh.Mesh(cat(3,XX,YY)) ;
        end
    % Initialize the levelset
        lvl = geo.Function(mesh.Nodes) ;
        lvlst = pkg.geometry.levelset.Mesh(mesh,lvl,'AutoBuildEdges',false) ;
    % Set as characteristic
        lvlst.Thickness = lvlst.defaultThickness ; % set chacteristic
    % Backup
        lvlst0 = lvlst ;
    % Build the velocity field
        V = (angle/nIt)*((mesh.Nodes-mean(mesh.Nodes,1))*[0 1;-1 0]) ;
    % Point tracers
        P = [] ; [mean(bbox,1)+[0 -.25].*range(bbox,1) 1] ;
    % Init display
        cla ; axis equal
        pl = plot(lvlst,'asmesh','EdgeColor','none') ;
        pt = patch('vertices',P,'Faces',(1:size(P,1))','marker','o','markeredgecolor','k') ;
        quiver(mesh.Nodes(:,1),mesh.Nodes(:,2),V(:,1),V(:,2)) ;
        drawnow ;
    % Updating loop
        profile on
        for it = 1:nIt
            if ~isempty(P)
                P(:,1:2) = P(:,1:2) + mesh.interpMat(P(:,1:2)+.0*(mesh.interpMat(P(:,1:2))*V))*V ; 
                pt.Vertices = P ;
            end
            lvlst = lvlst.advect(V) ;
            pl.Deformation = [0 0 1].*lvlst.LevelSetValues ;
            pl.CData = lvlst.LevelSetValues ;
            drawnow ;
        end
        profile off

%% COMPARE LEVEL-SETS
    cla ; axis equal
    %plot(lvlst0,'asmesh','FaceColor','none') ;
    %plot(lvlst,'asmesh','EdgeColor','none','FaceAlpha',.5) ;
    [~,edgP0] = lvlst0.buildEdges ;
    plot3(edgP0(:,1),edgP0(:,2),edgP0(:,1)*0+.5,'k') ;
    [~,edgP] = lvlst.buildEdges ;
    plot3(edgP(:,1),edgP(:,2),edgP(:,1)*0+.5,'r') ;



end
