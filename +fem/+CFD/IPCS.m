classdef IPCS < pkg.fem.AbstractFEM
% IPCS Transient CFD computation using the Incremental Pressure Correction Scheme
% see https://fenicsproject.org/pub/tutorial/sphinx1/._ftut1004.html#the-navier-stokes-equations
% - Navier stokes (velocity u, pressure p): 
%       rho*(du/dt + grad(u).u) = div(sigma(u,p)) + f
% - incompressible fluid: div(u) = 0
% - Non-Newtonian fluid: sigma(u,p) = 2*mu(d)*D(u) - p*Id
% - strain rate tensor D(u) = sym(grad(u)) ;
% - strain rate d = |D| ;
% - theta-scheme: u{n+theta} = theta*u* + (1-theta)*u{n} = theta*du* + u{n}
% - IPCS: three-step scheme (with N the outgoing boundary normal)
%       1) tentative velocity u*: v virtual velocity
%               int(rho*v.(u*-u{n})/dt*dx) 
%                   + int(rho*v.(u{n}.grad(u{n})*dx) 
%                   + int(D(v):sigma(u{n+theta},p{n})*dx)
%                   + int(v.N*p{n}*ds)
%                   - int(mu(d)*v.(grad(u{n+theta}).N)*ds)
%                   - int(v.f{n+1}*dx) = 0
%       2) pressure p{n+1}: q virtual pressure and dp = p{n+1} - p{n}
%               int(grad(q).grad(dp)*dx)*dt + int(q*div(u*)*dx) = 0
%       3) corrected velocity u{n+1}: du = u{n+1} - u*
%               int(v.du*dx) + dt*int(v.grad(dp)*dx) = 0
% - convention [u.grad(u)]_i = uj*dui_dxj

%% PROPERTIES
properties
% Master mesh
    Mesh pkg.geometry.mesh.Mesh 
% Function space meshes
    UMesh pkg.geometry.mesh.Mesh % velocity mesh
    PMesh pkg.geometry.mesh.Mesh % pressure mesh
% Nodal values
    u % velocity
    p % pressure
% Material properties
% - scalar: uniform distribution
% - vector: 
%   - [Mesh.nNodes 1] nodal values at the master mesh nodes
%   - [this.nGP 1] values at the gauss quadrature points
% - function handle @(d): function of the local strain rate d
    mu = 1 
    rho = 1
% Nodal forces: [UMesh.nNodes UMesh.nCoord]
    f = 0
% Boundary conditions
    Ubc pkg.fem.bc.Dirichlet
    Pbc pkg.fem.bc.Dirichlet
% Solver
    dt = []
    Theta = .5
    Solver = 'iterative'
    RelTol = 1e-6 
    MaxIt = 100 
end
properties %(Hidden)
% Constant objects
    InteriorQuadrature,BoundaryQuadrature,BoundaryNormals
    We,Wb
    Wu,Nu,Gu,gradU,divU,symGradU
    Wub,Nub,Gub,gradUb
    Np,Gp,gradP
    Npb,Gpb,gradPb
% Non-Fixed DOFs
    keepUdof
    keepPdof
% Material properties at qurature points
    MUi,MUb
    RHOi,RHOb
% FEM Matrices & preconditionners
    Add % strain rate
    Ads % boundary term
    A1,L1,U1,b1
    A2,L2,b2
    A3,L3,b3
end
    
%% METHODS
methods
    function this = IPCS(mesh,varargin)
    % Class constructor
        if nargin==0 ; error('A mesh must be provided') ; end
        this.Mesh = mesh ;
    % Build mesh function spaces
        this.UMesh = this.Mesh.setElementTypes(pkg.geometry.mesh.elements.LagrangeElement(this.Mesh.Elems.Types,2)) ;
        this.PMesh = this.Mesh.setElementTypes(pkg.geometry.mesh.elements.LagrangeElement(this.Mesh.Elems.Types,1)) ;
    % Initialize the results
        this.u = zeros(size(this.UMesh.Nodes)) ;
        this.p = zeros(this.PMesh.nNodes,1) ;
    % Process arguments
        if mod(numel(varargin),2) ; error('Wrong number of input arguments') ; end
        for arg = 1:2:numel(varargin)
            this.(varargin{arg}) = varargin{arg+1} ;
        end
    % Assemble the problem
        this.cst() ;
        this.lhs() ;
    end
    
    function dt = defaultTimeStep(this)
    % Return the default time step
        dt = .5*min(this.Mesh.elemSize(this.Mesh.Edges))/max(sqrt(sum(this.u.^2,2))) ;
    end
    
    function cst(this)
    % Compute the constant objects
        dim = this.UMesh.nCoord ;
    % Mesh integration
        intMesh = this.UMesh ;
        [ee,we,ie] = intMesh.integration() ; 
        this.InteriorQuadrature = [ie we ee] ;
        [eb,wb,ib,normal] = intMesh.boundaryIntegration() ; 
        this.BoundaryQuadrature = [ib wb eb] ;
        this.BoundaryNormals = normal ;
        this.We = diag(sparse(we)) ;
        this.Wb = diag(sparse(wb)) ;
    % Velocity interpolation
    % In the domain
        this.Wu = kron(speye(dim),this.We) ;
        this.Nu = kron(speye(dim),this.UMesh.interpMat(ee,ie)) ;
        this.Gu = this.UMesh.diffMat(ee,ie) ;
        this.gradU = kron(speye(dim),cat(1,this.Gu{:})) ; % gradU*u(:) = [du1_dx1 ; du1_dx2 ; du2_dx1 ; du2_dx2] ;
        this.divU = cat(2,this.Gu{:}) ;
    % On the boundary
        this.Wub = kron(speye(dim),this.Wb) ;
        this.Nub = kron(speye(dim),this.UMesh.interpMat(eb,ib)) ;
        this.Gub = this.UMesh.diffMat(eb,ib) ;
        this.gradUb = kron(speye(dim),cat(1,this.Gub{:})) ;
    % Symmetric part of the velocity gradient
    %   D(u) = sym(grad(u)) = grad(u) + grad(u)' ;
        this.symGradU = cellfun(@(mat)kron(speye(dim),mat),this.Gu(:)','UniformOutput',false) ;
        this.symGradU = .5*(this.gradU + cat(1,this.symGradU{:})) ;
    % Pressure interpolation
    % In the domain
        this.Np = this.PMesh.interpMat(ee,ie) ;
        this.Gp = this.PMesh.diffMat(ee,ie) ;
        this.gradP = cat(1,this.Gp{:}) ;
    % On the boundary
        this.Npb = this.PMesh.interpMat(eb,ib) ;
        this.Gpb = this.PMesh.diffMat(eb,ib) ;
        this.gradPb = cat(1,this.Gpb{:}) ;
    end
    
    function d = strainRate(this,ee,ie)
    % Return the strain rate evaluated at local coordinates
        if nargin<2 ; M = [this.gradU ; this.gradUb] ;
        else ; M = this.UMesh.gradMat(this.UMesh.nCoord,ee,ie) ;
        end
        G = M*this.u(:) ; % velocity gradient
        D = .5*(G + G(:,[1 3 2 4])) ; % symmetric part
        Ddev = D-[1/3 0 0 1/3].*sum(D.*[1 0 0 1],2) ; % Deviatoric part
        d = sqrt(sum(Ddev.^2,2)) ;
    end
    
    function [val,valBnd] = evalProp(this,prop)
    % Eval material properties at quadrature points
        if isscalar(prop) % uniform property
            val = prop ;
            valBnd = prop ;
        else % variable property
            if isa(prop,'function_handle') % function of the local strain rate
                val = reshape(prop(this.strainRate()),[],1) ;
            elseif numel(prop)==this.PMesh.nNodes % nodal values given on the pressure mesh
                val = [this.Np ; this.Npb]*prop(:) ;
            elseif numel(prop)==this.UMesh.nNodes % nodal values given on the velocity mesh
                val = [this.Nu(1:end/2,1:end/2) ; this.Nub(1:end/2,1:end/2)]*prop(:) ;
            elseif numel(prop)==numel(ie) % given at quadrature points
                val = prop(:) ;
            else
                error('Wrong format for the propertie values') ;
            end
        % Split interior and boundary values
            nQi = size(this.InteriorQuadrature,1) ;
            valBnd = val(1:nQi) ;
            val = val(nQi+1:end) ;
        end
    end

    function lhs(this,onlyA1)
    % Assemble the Left-Hand-Side objects
    % May Compute only A1 (when rho and/or mu has changed)
        if nargin<2 ; onlyA1 = false ; end
    % Boundary conditions
        if ~onlyA1
        % Velocity
            this.keepUdof = true(size(this.u)) ;
            for bc = this.Ubc(:)'
                [u0,fixDOF] = bc.apply(this.UMesh) ;
                this.keepUdof = this.keepUdof & ~fixDOF ;
                this.u(fixDOF,:) = u0(fixDOF,:) ;
            end
        % Pressure
            this.keepPdof = true(size(this.p)) ;
            for bc = this.Pbc(:)'
                [p0,fixDOF] = bc.apply(this.PMesh) ;
                this.keepPdof = this.keepPdof & ~fixDOF ;
                this.p(fixDOF) = p0(fixDOF) ;
            end
        end
    % Material properties
        [this.RHOi,this.RHOb] = this.evalProp(this.rho) ;
        [this.MUi,this.MUb] = this.evalProp(this.mu) ;
    % 1) tentative velocity u*: v virtual velocity and u{n+1/2} = (1-theta)*u{n} + theta*u* 
        dim = this.UMesh.nCoord ;
        this.A1 = sparse(numel(this.UMesh.Nodes),numel(this.UMesh.Nodes)) ;
        % int(rho*v.(u*-u{n})/dt*dx) 
            if ~isempty(this.dt)
                M = this.Nu'*(this.Wu*diag(sparse(this.RHOi/this.dt)))*this.Nu ;
                this.A1 = this.A1 + M ;
            end
        % + int(D(v):sigma(u{n+1/2},p{n})*dx)
        %   with sigma(u,p) = 2*mu(d)*D(u) - p*Id
            this.Add = this.symGradU'*kron(speye(dim),this.Wu*diag(sparse(2*this.MUi)))*this.symGradU ;
            this.A1 = this.A1 + this.Theta*this.Add ;
        % - int(mu(d)*v.(grad(u{n+1/2}).N)*ds) % [grad(u).n]_i = nj*duj_dxi
            nn = cellfun(@(ni)kron(speye(dim),diag(sparse(ni))),num2cell(this.BoundaryNormals,1),'UniformOutput',false) ;
            nn = cat(2,nn{:}) ;
            this.Ads = this.Nub'*(this.Wub*diag(sparse(this.MUb)))*(nn*this.gradUb) ;
            this.A1 = this.A1 - this.Theta*this.Ads ;
        % Apply BCs
            this.A1 = this.A1(this.keepUdof,this.keepUdof) ;
        % Preconditionning
            [this.L1,this.U1] = ilu(this.A1) ;
    % Compute only A1 ? (when rho and/or mu has changed)
        if onlyA1 ; return ; end
    % 2) pressure p{n+1}: q virtual pressure and dp = p{n+1} - p{n}
        %   int(grad(q).grad(dp)*dx)*dt + int(q*div(u*)*dx) = 0
            this.A2 = this.gradP'*this.Wu*this.gradP ;
            this.A2 = this.A2(this.keepPdof,this.keepPdof) ;
            this.L2 = ichol(this.A2) ;
    % 3) corrected velocity u{n+1}: du = u{n+1} - u*        
        %   int(v.du*dx) + dt*int(v.grad(dp)*dx) = 0
            this.A3 = this.Nu'*this.Wu*this.Nu ; 
            this.A3 = this.A3(this.keepUdof,this.keepUdof) ;
            this.L3 = ichol(this.A3) ;
    end

    function steady(this)
    % Steady-state solution
    end

    function step(this)
    % Stepping scheme
        dim = this.UMesh.nCoord ;
    % Set default time step
        if isempty(this.dt) ; this.dt = this.defaultTimeStep ; end
    % 1) tentative velocity u*: v virtual velocity and u{n+1/2} = (1-theta)*u{n} + theta*u* 
        % int(rho*v.(u*-u{n})/dt*dx) 
            this.b1 = zeros(numel(this.u),1) ; 
        % + int(rho*v.(u{n}.grad(u{n})*dx) 
            uGradU = sum((this.Nu*this.u(:)).*reshape(this.gradU*this.u(:),[],dim),2) ; % [u.grad(u)]_i = uj*dui_dxj
            this.b1 = this.b1 - this.Nu'*(this.Wu*diag(sparse(this.RHOi)))*uGradU(:) ;
        % + int(D(v):sigma(u{n+1/2},p{n})*dx)
        %   with sigma(u,p) = 2*mu(d)*D(u) - p*Id
            this.b1 = this.b1 - this.Add*this.u(:) ;
            this.b1 = this.b1 + this.symGradU'*kron(speye(dim),this.Wu)*kron(reshape(speye(dim),[],1),this.Np*this.p(:)) ;
        % + int(v.N*p{n}*ds)
            this.b1 = this.b1 - this.Nub'*diag(sparse(this.BoundaryNormals(:)))*repmat(this.Wb,dim,1)*(this.Npb*this.p(:)) ;
        % - int(mu(d)*v.(grad(u{n+1/2}).N)*ds) % [grad(u).n]_i = nj*duj_dxi
            this.b1 = this.b1 + this.Ads*this.u(:) ;
        % - int(v.f{n+1}*dx)
            if any(this.f~=0)
                this.b1 = this.b1 + this.Nu'*(this.Wu*(this.Nu*this.f(:))) ;
            end
        % Solve
            this.b1 = this.b1(this.keepUdof) ;
            switch this.Solver
                case 'direct'
                    du = this.A1\this.b1 ;
                case 'iterative'
                    [du,flag] = gmres(this.A1,this.b1,[],this.RelTol,this.MaxIt,this.L1,this.U1) ;
                    %[u(keepUdof),flag] = bicgstabl(A1,b1,[],100,L1,U1,u(keepUdof)) ;
                    if flag>0 ; warning('GMRES failed') ; end
            end
            this.u(this.keepUdof) = this.u(this.keepUdof) + du ;
    % 2) pressure p{n+1}: q virtual pressure and dp = p{n+1} - p{n}
        %   int(grad(q).grad(dp)*dx)*dt + int(q*div(u*)*dx) = 0
            this.b2 = - (1/this.dt)*this.Np'*this.We*(this.divU*this.u(:)) ;
        % Solve
            this.b2 = this.b2(this.keepPdof) ;
            dp = zeros(size(this.p)) ;
            switch this.Solver
                case 'direct'
                    dp(this.keepPdof) = this.A2\this.b2 ;
                case 'iterative'
                    [dp(this.keepPdof),flag] = pcg(this.A2,this.b2,this.RelTol,this.MaxIt,this.L2,this.L2') ;
                    %[dp(keepPdof),flag] = bicgstabl(A2,b2,reltol,maxIt,L2,L2') ;
                    if flag>0 ; warning('PCG failed') ; end
            end
            this.p = this.p + dp ;
   % 3) corrected velocity u{n+1}: du = u{n+1} - u*
        %   int(v.du*dx) + dt*int(v.grad(dp)*dx) = 0 
            this.b3 = -this.dt*(this.Nu'*(this.Wu*(this.gradP*dp))) ;
        % Solve
            this.b3 = this.b3(this.keepUdof) ;
            switch this.Solver
                case 'direct'
                    du = this.A3\this.b3 ;
                case 'iterative'
                    [du,flag] = pcg(this.A3,this.b3,this.RelTol,this.MaxIt,this.L3,this.L3') ;
                    if flag>0 ; warning('PCG failed') ; end
            end
            this.u(this.keepUdof) = this.u(this.keepUdof) + du ;
    end
end


%% EXAMPLES
methods (Static)
    function channelFlow
    %% FLOW IN A CHANNEL
        clc,clearvars
        L = [1 .3] ; dx = min(L)/10 ;
        mu = 1e2 ; rho = 1 ; 
        dP = 10 ; dt = dx/10 ;
        mesh = pkg.geometry.levelset.Rectangle(L.*[0;1]).mesh(dx) ;
        profile on
        Ubc = [pkg.fem.bc.Dirichlet([NaN 0],[0 0]) pkg.fem.bc.Dirichlet([NaN L(2)],[0 0])] ;
        Pbc = [pkg.fem.bc.Dirichlet([0 NaN],dP) pkg.fem.bc.Dirichlet([L(1) NaN],0)] ;
        FEM = pkg.fem.CFD.IPCS(mesh,'rho',rho,'mu',mu,'Ubc',Ubc,'Pbc',Pbc,'dt',dt) ;
        clf ; axis equal ; pl = plot(FEM.PMesh) ; q = quiver(FEM.UMesh.Nodes(:,1),FEM.UMesh.Nodes(:,2),FEM.u(:,1),FEM.u(:,2)) ;
        T=0 ;
        while T<1000*dt
            tic ; FEM.step ; toc
            T = T+FEM.dt ;
            pl.CData = FEM.p ;
            q.UData = FEM.u(:,1) ; q.VData = FEM.u(:,2) ;
            drawnow ;
        end
        profile off
        
    end
end

end

%% UNITARY TESTS
function tests
end

