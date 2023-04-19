%% COMPUTE PLANE WAVE PROPAGATION IN WAVEGUIDES

%% WAVES IN AN HOMOGENEOUS PRISMATIC WAVEGUIDE
% Parameters
    % Geometry
        L = [4 .5] ; % unit cell size (mm)
        dx = .24 ; % element size (mm)
        perVec = diag(L) ; % periodicity vectors
    % Material
        E = .5e3*(1+.00001i) ; % Young modulus (MPa)
        G = E/2.6 ; % shear modulus (MPa)
        nu = E./(2*G)-1 ;
        rho = 1200e-12 ; % material density (tons/mm^3)
    % Wave 
        dir = [NaN NaN 1] ; % wave propagation direction [Nan: no periodicity]
        freq = logspace(3,6,300) ; % wave frequency
        nModes = 50 ; % number of extracted wave modes
    % Display
        plotType = 'k' ; % 'c', 'k' , 'l' or 'g'
        logScale = true ;
        gammaMax = .1 ; inf ;
% Build the mesh
    mesh = pkg.geometry.mesh.GridMesh(L,dx) ;
    %if numel(L)==2 ; mesh.setElementTypes(pkg.geometry.mesh.elements.Quad8) ; end
% Display the mesh
    clf ; axis equal tight ; view([30 30]) ;
    plot(mesh) ;
% Build the FEM matrices
    C = stiffness(E,G) ;
    [K00,K0i,Kij,M,P] = FEM(mesh,C,rho,[],perVec) ;
% Compute the wavenumbers and modes
    [K,U,omega] = bloch(mesh,K00,K0i,Kij,M,P,freq,dir,nModes) ;
% DISPLAY RESULTS
    clf ;
% Theoretical wavenumbers
        b = L(1) ; h = L(2) ;
    % Longitudinal wave
        S = b*h ;
        kl = omega.*sqrt(rho/E) ;
        plot3(freq,real(kl),imag(kl),'-','linewidth',1) ;
    % Bending waves
        I11 = b^3*h/12 ;
        I22 = b*h^3/12 ;
        kb1 = sqrt(omega).*(rho*S/E/I11)^.25 ;
        kb2 = sqrt(omega).*(rho*S/E/I22)^.25 ;
        plot3(freq,real(kb1),imag(kb1),'-','linewidth',1) ;
        plot3(freq,real(kb2),imag(kb2),'-','linewidth',1) ;
    % ... with transverse shear correction
        xi = pi^2/12 ; p = rho/E ; q = rho/xi/G ;
        kb1s1 = omega.*sqrt((p+q) - sqrt((p-q)^2 + 4./(omega.^2).*(rho*S/E/I11)))/sqrt(2) ;
        plot3(freq,real(kb1s1),imag(kb1s1),'-','linewidth',1) ;
        kb1s2 = omega.*sqrt((p+q) + sqrt((p-q)^2 + 4./(omega.^2).*(rho*S/E/I11)))/sqrt(2) ;
        plot3(freq,real(kb1s2),imag(kb1s2),'-','linewidth',1) ;
        kb2s1 = omega.*sqrt((p+q) - sqrt((p-q)^2 + 4./(omega.^2).*(rho*S/E/I22)))/sqrt(2) ;
        plot3(freq,real(kb2s1),imag(kb2s1),'-','linewidth',1) ;
        kb2s2 = omega.*sqrt((p+q) + sqrt((p-q)^2 + 4./(omega.^2).*(rho*S/E/I22)))/sqrt(2) ;
        plot3(freq,real(kb2s2),imag(kb2s2),'-','linewidth',1) ;
    % Torsion wave
        a = interp1([1 1.5 2 2.5 3 4 5 6 10 1000],...
                    [.141 .196 .229 .249 .263 .281 .291 .299 .312 .333],...
                    max(b,h)/min(b,h)) ; % torsion correction, see en.wikipedia.org/wiki/Torsion_constant
        It = a*max(b,h)*min(b,h)^3 ;
        J = I11+I22 ;
        kt = omega.*(rho*J/G/It)^.5 ;
        plot3(freq,real(kt),imag(kt),'-','linewidth',1) ;
% Computed wavenumbers
    FREQ = freq + 0*K ;
    Kdir = K(:).*dir ;
    plothandle = plot3(FREQ(:),real(K(:)),imag(K(:)),'+k','markersize',4,'linewidth',.1) ;
% Set plot options
    setPlot(gammaMax,plotType,logScale) ;
    waveModeAnimation(mesh,Kdir,U,plothandle)
    

%% WAVES IN A TWO-PHASE INFINITE 3D SOLID
% Parameters
    % Geometry
        dx = .1 ; % voxel size (mm)
        L = 5*dx*[1 1 1]*(1-eps) ; % unit cell size (mm)
        subdiv = 2 ; % voxel subdivision
        perVec = diag(L) ; % periodicity vectors
    % Random phase distribution
        volfrac = .15 ;
        n = prod(ceil(L/dx)) ;
        softcell = randsample(n,floor(n*volfrac)) ; [] ; % indices
        issoft = full(sparse(softcell,1,true,n,1)) ; % logical array
    % Voxel subdivision
        issoft = reshape(issoft,ceil(L/dx)) ;
        if subdiv>1
            issoft = repelem(issoft,subdiv,subdiv,subdiv) ;
            dx = dx/subdiv ;
        end
    % Material
        E = 4e3*(1+.001i)*~issoft + 1e3*(1+.001i)*issoft ; % Young modulus (MPa)
        G = E/2.6 ; % shear modulus (MPa)
        rho = 1200e-12*~issoft + 1100e-12*issoft ; % material density (tons/mm^3)
    % Wave 
        dir = [1 0 0] ; % wave propagation direction [Nan: no periodicity]
        freq = logspace(5,7,100) ; linspace(.1e6,10e6,20) ; logspace(5,7,20) ; % wave frequency
        nModes = 10 ; % number of extracted wave modes
    % Display
        plotType = 'c' ; % 'c', 'k' , 'l' or 'g'
        logScale = true ;
        gammaMax = .1 ; inf ;
% Build the mesh
    mesh = pkg.geometry.mesh.GridMesh(L,dx) ;
% Display the mesh
    soft = copy(mesh) ; soft.Elems = mesh.Elems.subpart(issoft) ;
    clf ; axis equal tight ; view([30 30]) ;
    plot(mesh,'VisibleEdges','none','FaceAlpha',0) ;
    plot(soft,'FaceColor',[1 1 1]*.5,'FaceAlpha',0.5,'HighlightBoundaryEdges',false) ;
% Build the FEM matrices
    C = stiffness(E,G) ;
    [K00,K0i,Kij,M,P] = FEM(mesh,C,rho,[],perVec) ;
% Compute the wavenumbers and modes
    [K,U,omega] = bloch(mesh,K00,K0i,Kij,M,P,freq,dir,nModes) ;
% DISPLAY RESULTS
    clf ;
% Theoretical wavenumbers
    rhom = mean(rho(:)) ;
    Cvoigt = mean(C,3) ;
    Creuss = 1./mean(1./C,3) ;
    % Longitudinal wave
        klr = omega.*sqrt(rhom/Creuss(1,1)) ;
        klv = omega.*sqrt(rhom/Cvoigt(1,1)) ;
        patch([freq flip(freq)],[real(klr) flip(real(klv))],[imag(klr) flip(imag(klv))],'r')
    % Transverse wave
        ktr = omega.*sqrt(rhom/Creuss(6,6)) ;
        ktv = omega.*sqrt(rhom/Cvoigt(6,6)) ;
        patch([freq flip(freq)],[real(ktr) flip(real(ktv))],[imag(ktr) flip(imag(ktv))],'b')
    set(findobj(gca,'type','patch'),'edgecolor','none','facealpha',.25)
% Computed wavenumbers
    FREQ = freq + 0*K ;
    Kdir = K(:).*dir ;
    plothandle = plot3(FREQ(:),real(K(:)),imag(K(:)),'+k','markersize',4,'linewidth',.1) ;
    %plot3(freq,real(K),imag(K),'+k','markersize',4,'linewidth',.1) ;
% Set plot options
    setPlot(gammaMax,plotType,logScale) ;
    waveModeAnimation(mesh,Kdir,U,plothandle) ;


%% MATERIAL STIFFNESS MATRIX
function C = stiffness(E,G)
    E = reshape(E,1,1,[]) ; G = reshape(G,1,1,[]) ; O = zeros(size(E)) ;
    nu = E./(2*G)-1 ;
    C = E./(1+nu)./(1-2*nu).*[...
                              1-nu nu nu O O O ; ...
                              nu 1-nu nu O O O ; ...
                              nu nu 1-nu O O O ; ...
                              O O O .5-nu O O ; ...
                              O O O O .5-nu O ; ...
                              O O O O O .5-nu ; ...
                             ] ;
end

%% SYSTEM MATRICES
% u(x,t) = U(x)*exp(1i*omega*t - 1i*k*x)
% (K00 + k.K1 + k^2.K2 - w^2.M).U = 0

function [K00,K0i,Kij,M,P] = FEM(mesh,C,rho,sig,perVec)

    % Mesh interpolation
    [ee,we,ie] = mesh.integration() ;
    nQP = numel(we) ;
    W = diag(sparse(we)) ;
    O = sparse(nQP,mesh.nNodes) ;
    N = mesh.interpMat(ee,ie) ;
    D = mesh.diffMat(ee,ie) ; 
    [D{end+1:3}] = deal(O) ;

    % EPS = [E11;E22;E33;2E13;2E23;2E12] = B*U = (B0 + sum(B{i}*k_i))*u = [
    B0 = [D{1} O O ; O D{2} O ; O O D{3} ; D{3} O D{1} ; O D{3} D{2} ; D{2} D{1} O] ;
    Bi = cell(3,1) ;
    Bi{1} = -1i*[N O O ; O O O ; O O O ; O O N ; O O O ; O N O] ;
    Bi{2} = -1i*[O O O ; O N O ; O O O ; O O O ; O O N ; N O O] ;
    Bi{3} = -1i*[O O O ; O O O ; O O N ; N O O ; O N O ; O O O] ;

    % DISTRIBUTE MATERIAL STIFFNESS & WEIGHTS
    if size(C,3)==1 % homogeneous solid
        CW = kron(sparse(C),W) ;
    else % heterogeneous
        Cqp = C(:,:,ie) ; % material stiffness at quadrature points
        Cqp = Cqp.*reshape(we,1,1,[]) ; % times weights
        iii = (0:5)'*nQP + zeros(1,6) + reshape(1:nQP,1,1,[]) ;
        jjj = zeros(6,1) + (0:5)*nQP + reshape(1:nQP,1,1,[]) ;
        CW = sparse(iii(:),jjj(:),Cqp(:),6*nQP,6*nQP) ;
    end

    % STIFFNESS MATRICES
    K00 = B0'*CW*B0 ; % zero-order matrix
    K0i =  cell(3,1) ; % first-order matrices
    Kij = cell(3,3) ; % second order matrices
    for ii = 1:3 
        K0i{ii} = B0'*CW*Bi{ii} + Bi{ii}'*CW*B0 ;
        for jj = 1:3
            Kij{ii,jj} = Bi{ii}'*CW*Bi{jj} ; % Kij{ii,jj} = Kij{jj,ii}'
        end
    end

    % MASS MATRICES
    if numel(rho)==1
        rhoW = rho*W ;
    else
        rhogp = rho(ie) ;
        rhoW = diag(sparse(rhogp(:).*we)) ;
    end
    M = (N'*rhoW*N) ;
    M = blkdiag(M,M,M) ;
    
    % PRESTRESS INFLUENCE
    if any(abs(sig)>eps)
    % Stress components on quadrature points
        if size(sig,1)==1 % homogeneous stress
            sigW = we.*sig ;
        elseif size(sig,1)==mesh.nNodes % nodal stress
            sigW = W*(N*sig) ;
        else % element stress
            sigW = W*sig(ie,:) ;
        end
    % Distibute as a sparse matrix
        sigW = permute(sigW(:,[1 6 4;6 2 5;4 5 3]),[2 3 1]) ; % (3,3,nQP) 
        iii = (0:2)'*nQP + zeros(1,3) + reshape(1:nQP,1,1,[]) ;
        jjj = zeros(3,1) + (0:2)*nQP + reshape(1:nQP,1,1,[]) ;
        sigW = sparse(iii(:),jjj(:),sigW(:),3*nQP,3*nQP) ;
    % Displacement gradients [u1,1 ; u1,2 ; u1,3 ; ... ; u3,3]
        G0 = cat(1,D{:}) ;
        Gk = -1i*N ;
    % Update stiffness matrices
        k00 = G0'*sigW*G0 ;
        K00 = K00 + blkdiag(k00,k00,k00) ;
        for ii = 1:3 
            Gi = kron(sparse(ii,1,1,3,1),Gk) ;
            k0i = G0'*sigW*Gi ;
            k0i = k0i + k0i' ; % symmetric hermitian
            K0i{ii} = K0i{ii} + blkdiag(k0i,k0i,k0i) ;
            for jj = 1:3
                Gj = kron(sparse(jj,1,1,3,1),Gk) ;
                kij = Gi'*sigW*Gj ;
                kij = kij + kij' ; % symmetric hermitian
                Kij{ii,jj} = Kij{ii,jj} + blkdiag(kij,kij,kij) ;
            end
        end
    end

    % PERIODICITY MATRICES
    P = cell(size(perVec,1),1) ;
    if ~isempty(perVec)
        [P{:}] = mesh.perNodeMat(perVec) ;
    end
    [P{end+1:3}] = deal(speye(mesh.nNodes)) ; 
    
end


%% BLOCH WAVE COMPUTATION
% (K00 + mu.K1 + mu^2.K2 - w^2.M).U = 0
% V = mu.U
% ([w^2.M-K0 0;0 K2] - mu.[K1 K2;K2 0])[U;V] = [0;0]

function [K,U,omega] = bloch(mesh,K00,K0i,Kij,M,P,freq,dir,nModes)

    omega = 2*pi*freq ;
    nDir = size(dir,1) ;
    nFreq = numel(omega) ;
    nDOFs = size(M,1) ;
    nDims = nDOFs/mesh.nNodes ;

    vK0i = cellfun(@(M)reshape(M,[],1),K0i,'uni',false) ;
    vK0i = cat(2,vK0i{:}) ;
    vKij = cellfun(@(M)reshape(M,[],1),Kij,'uni',false) ;
    vKij = cat(2,vKij{:}) ;

    K = NaN(nModes,nFreq,nDir) ;
    U = NaN(mesh.nNodes,nDims,nModes,nFreq,nDir) ;
    wtbr = waitbar(0,'Bloch Wave Computation..') ;
    for dd = 1:nDir
    % Full wave direction vector
        n = dir(dd,:) ;
        isWaveDir = ~isnan(n) ;
        n(~isWaveDir) = 0 ; % avoid multiplication by NaNs
    % Periodocity
        T = speye(mesh.nNodes) ;
        for dim = 1:nDims
            if isWaveDir(dim) 
                T = P{dim}*T ; 
            end
        end
        T(:,sum(T,1)==0) = [] ; % delete unused DOFs
        T = kron(speye(nDims),T) ;
        T2 = blkdiag(T,T) ; % for the polynomial EVP
    % Projected stiffnesses
        nn = n(:).*n(:)' ; 
        K1 = reshape(vK0i*n(:),size(M)) ;
        K2 = reshape(vKij*nn(:),size(M)) ;
        I = speye(size(M)) ; % K2 ;
    % For each frequency..
        for ww = 1:nFreq
            K0 = omega(ww)^2*M-K00 ;
            A = blkdiag(K0,I) ; 
            B = [K1 K2;I sparse(nDOFs,nDOFs)] ;
            [uu,ku] = eigs(T2'*A*T2,T2'*B*T2,nModes,'sm') ;
            K(:,ww,dd) = diag(ku) ;
            U(:,:,:,ww,dd) = reshape(T*uu(1:size(T,2),:),mesh.nNodes,nDims,nModes) ;
            wtbr = waitbar((ww+nFreq*(dd-1))/nFreq/nDir,wtbr) ;
        end
    end
    delete(wtbr) ;

end

%% Display
function setPlot(gammaMax,plotType,logScale)

    % Get lines
    lines = findobj(gca,'type','line')' ;
    lines = [lines findobj(gca,'type','patch')'] ;

    % CullgGamma>gammaMax
    for ll = lines 
        ki = ll.YData + 1i*ll.ZData ;
        valid = ones(size(ki)) ; 
        valid(abs(imag(ki)./real(ki))>gammaMax) = NaN ; 
        ll.YData = ll.YData.*valid ;
        ll.ZData = ll.ZData.*valid ;
    end

    % Wavenumber, phase velocity, etc
    xlabel 'Frequency (Hz)'
    switch plotType
        case 'c'
            ylabel 'Phase Velocity (m/s)'
            for ll = lines 
                omega = 2*pi*ll.XData ;
                ki = ll.YData + 1i*ll.ZData ;
                ll.YData = 1e-3*real(omega./ki) ; 
                ll.ZData = 1e-3*imag(omega./ki) ; 
            end
        case 'k'
            ylabel 'Wavenumber (rad/mm)'
        case 'g'
            ylabel 'Wavenumber (rad/mm)'
            zlabel 'Spatial decay'
            for ll = lines 
                ki = ll.YData + 1i*ll.ZData ;
                ll.ZData = imag(ki)./real(ki) ; 
            end
        case 'l'
            ylabel 'Wavelength (mm)'
            for ll = lines 
                ki = ll.YData + 1i*ll.ZData ;
                ll.YData = real(2*pi./ki) ; 
                ll.ZData = imag(2*pi./ki) ; 
            end
    end

    % Scale change
    if logScale
        for ll = lines 
            ll.YData = abs(ll.YData) ; 
            ll.ZData = abs(ll.ZData) ; 
        end
        set(gca,'xscale','log','yscale','log','zscale','log')
    end

end

function waveModeAnimation(mesh,Kdir,U,plothandle)
    tag = 'waveAnim' ;
    amp = 10/100./norm(range(mesh.Nodes,1)) ;
    timerPeriod = .05 ;
    animFreq = .5  ;
    Le = 3.5*norm(range(mesh.boundingBox,1)) ;
    de = median(mesh.elemSize(mesh.Edges)) ;

    % Delete all previous objects
    delete(timerfindall('tag',tag)) ;
    delete(findall(0,'tag',tag)) ;
    
    % Create the cursor
    fig = gcf ; ax = gca ;
    curs = patch('vertices',NaN(1,3)...
                ,'faces',[1 1 1]...
                ,'marker','.'...
                ,'EdgeColor','r'...
                ,'MarkerSize',20 ...
                ,'tag',tag) ;
    
    % Data points
    pts = [plothandle.XData(:) plothandle.YData(:) plothandle.ZData(:)] ;
                            
    % Animation figure
    figure('tag',tag) ;
    axis equal tight off ; view([30 30]) ;
    
    % Display a 3D mesh
    mesh3D = mesh ;
    for cc = mesh.nCoord+1:3
        mesh3D = mesh3D.extrude(Le*full(sparse(1,mesh3D.nCoord+1,1,1,3)),ceil(Le/de)) ;
    end
    pl = plot(mesh3D) ;
    
    % 3D mode
    mode3D = @(pp)repmat(U(:,:,pp),[mesh3D.nNodes/mesh.nNodes 1]).*exp(-1i*sum(Kdir(pp,:).*mesh3D.Nodes,2,'omitnan')) ;
    curs.UserData = mode3D(1) ;
    
    % Cursor function
    fig.WindowButtonMotionFcn = @(src,evt)cellfun(@(c)c(src,evt),{...
                    @(src,evt)set(curs,'vertices',pts(closestToMouse(pts,ax),:)) ...
                    @(src,evt)set(curs,'UserData',mode3D(closestToMouse(pts,ax))) ...
                                },'uni',false);

    % Timer that makes the wave phase turn                        
    startTime = tic ;
    defShape = @()real(curs.UserData*exp(2i*pi*toc(startTime)*animFreq)) ;
    timerfcn = @(src,evt)set(pl,'Deformation',defShape(),'CData',sqrt(sum(defShape().^2,2))) ;
    ti = timer('Period',timerPeriod,'ExecutionMode','FixedRate','TimerFcn',timerfcn,'tag',tag) ;
    start(ti)

    % Delete timer on figure close
    addlistener(pl.Parent,'ObjectBeingDestroyed',@(src,evt)delete(ti)) ;
    addlistener(pl.Parent,'ObjectBeingDestroyed',@(src,evt)set(fig,'WindowButtonMotionFcn',[])) ;

end

    
function [p,d] = closestToMouse(pts,ax)
% Get the closest point to the mouse pointer
    if nargin<2 ; ax = gca ; end
% Mouse position in axes 
    line = ax.CurrentPoint ;
% Set to limits if rotation perpendicular
    lims = [ax.XLim(:) ax.YLim(:) ax.ZLim(:)] ;
    line(isinf(line)) = lims(isinf(line)) ;
% Set scales if needed
    islog = ismember({ax.XScale ax.YScale ax.ZScale},'log') ;
    line(:,islog) = log10(line(:,islog)) ;
    lims(:,islog) = log10(lims(:,islog)) ;
    pts(:,islog,:) = log10(pts(:,islog,:)) ;
% Normalize with limits
    m = min(lims,[],1) ;
    R = range(lims,1) ;
    line = (line-m)./R ;
    pts = (pts-m)./R ;
% Distance from line to all points
    u = diff(line,1,1)./sqrt(sum(diff(line,1,1).^2,2)) ;
    v = pts-line(1,:) ;
    t = sum(v.*u,2) ;
    d = sqrt(abs(sum(v.^2,2) - t.^2)) ;
    [d,p] = min(d(:))  ;
end










