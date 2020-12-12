%% SIMULTANEOUS SPACE-TIME DISCRETIZATION OF HEAT TRANSFER EQUATIONS

%% Problem Definition: 1D bar, moving point heat source

% Material Parameters
    L = 1 ; % m bar length
    S = .8e-3 ; % m² bar section
    Cp = 50 ; % Jouls/Kelvins/kg specific heat capacity
    lmbda = 1600 ; % Watts/Kelvins thermal conductivity
    rho = 800 ; % kg/m^3 density
    Tlim = [0 10] ; 'auto' ; % for caxis
% Source parameters
    Qs = 250 ; % Watts source heat flux
    T = 1 ; % sec time period
    dts = 0.001 ; % sec time increment
    ts = (0:dts:T*0.5)' ; % time samples
    ys = ... L/3*(1-ts/T) + 2*L/3*ts/T ... linear
         ... ts*0 + L/2 ... fixed
         ... L/2*(1+0.5*sign(sin(ts/T*2*pi*20))) ... alternate
          L/3*(1 + mod(4*ts/T,1)) ... saw
         ... L*0.5*(1+0.5*sin(ts*2*pi*5/T)) ... sinusoidal
         ; 
    Xs = [ts(:) ys(:)] ; % [t y] source position

    
%% Usual problem solving
% int_X[rho.Cp.T*.dT/dt] + int_X[dT*/dx.lambda.dT/dx] = int_X[T*.Qs]
% Parameters
    y = linspace(0,L,1000)' ;
    t = linspace(0,T,1000)' ;
    scheme = 'implicit' ; % 'explicit' or 'implicit'
    display = false ;
% Create the mesh
    mesh_1D = pkg.geometry.mesh.Mesh(y) ;
% Problem matrices
    [E,w,ie] = mesh_1D.integration ; W = spdiags(w(:),0,numel(w),numel(w)) ;
    N = mesh_1D.interpMat(E,ie) ;
    G = mesh_1D.gradMat(E,ie) ; G = G{1} ;
% Assemble
    K = S*lmbda*G'*W*G ; K = 0.5*(K+K') ;
    M = S*Cp*rho*N'*W*N ; M = 0.5*(M+M') ;
% Display
    if display
        clf
        pl = plot(mesh_1D) ;
    end
% Source path interpolation
    ys = interp1(Xs(:,1),Xs(:,2),t,'linear',Inf) ;
    Ns = mesh_1D.interpMat(ys) ;
    Ns(ys>L | ys<0,:) = 0 ;
% Loop over time range
    temp = zeros(mesh_1D.nNodes,numel(t)) ;
    dt = mean(diff(t)) ;
    Mt = M/dt ;
    wtbr = waitbar(0,'Transient Analysis..') ;
    wtbrTime = tic ;
    for tt = 2:numel(t)
    % Second-hand side vector
        f = Qs*Ns(tt,:)' ;
    % Updating
        switch scheme
            case 'explicit'
                temp(:,tt) = temp(:,tt-1) + (M\(f-K*temp(:,tt-1)))*dt ;
            case 'implicit'
                temp(:,tt) = (K+Mt)\(f + Mt*temp(:,tt-1)) ;
        end
    % Update display
        if toc(wtbrTime)>0.1
            wtbr = waitbar(tt/numel(t),wtbr) ; 
            wtbrTime = tic ; 
        end
        if display
            pl.NodeCoordinates = [y temp(:,tt)] ;
            drawnow ;
        end
    end
    delete(wtbr) ;
% Display the final result
    clf reset
    [TT,YY] = meshgrid(t,y) ;
    surf(TT,YY,temp) ;
    shading interp
    %title(['Space only, ' scheme])
    xlabel 'Time (sec)'
    ylabel 'Position (m)'
    colorbar
    axis equal,axis tight
    colormap(jet(16))
    caxis(Tlim)
   
    
%% SPACE-TIME DISCRETIZATION, Uniform Mesh
% int_(X,t)[rho.Cp.T*.dT/dt] + int_(X,t)[dT*/dx.lambda.dT/dx] = int_(X,t)[T*.Qs]
% Parameters
    y = linspace(0,L,100)' ;
    t = linspace(0,T,100)' ;
% Build the mesh
    [TT,YY] = meshgrid(t,y) ;
    mesh_2D = pkg.geometry.mesh.Mesh(cat(3,TT,YY)) ;
    %mesh_2D.setElementTypes(pkg.geometry.mesh.elements.Quad8) ;
% Problem matrices
    [E,w,ie] = mesh_2D.integration ; W = spdiags(w(:),0,numel(w),numel(w)) ;
    N = mesh_2D.interpMat(E,ie) ;
    G = mesh_2D.gradMat(E,ie) ; Gy = G{2} ; Gt = G{1} ;
% Assemble (DO NOT SYMMETRIZE)
    K = S*lmbda*Gy'*W*Gy + S*Cp*rho*N'*W*Gt ;
% Second-hand term
    dt = [ts(2)-ts(1) ; (ts(3:end)-ts(1:end-2))/2 ; ts(end)-ts(end-1)] ;
    f = Qs*mesh_2D.interpMat(Xs)'*dt ;
% Boundary conditions
    keepDOF = mesh_2D.Nodes(:,1)>eps ;
    K = K(keepDOF,keepDOF) ;
    f = f(keepDOF) ;
% Solve
    temp_st = zeros(mesh_2D.nNodes,1) ;
    temp_st(keepDOF) = K\f ;
% Display
    clf reset
    pl = plot(mesh_2D) ;
    pl.FaceColor = 'interp' ;
    pl.Faces.FaceVertexCData = temp_st ;
    pl.VisibleEdges = 'none' ;
    pl.NodeCoordinates = mesh_2D.Nodes.*[1 1] ;
    %title 'Space-Time, Uniform mesh'
    xlabel 'Time (sec)'
    ylabel 'Position (m)'
    colorbar
    axis equal,axis tight
    colormap(jet(16))
    caxis(Tlim)
   
    
%% SPACE-TIME DISCRETIZATION, Unstructured Mesh
% int_(X,t)[rho.Cp.T*.dT/dt] + int_(X,t)[dT*/dx.lambda.dT/dx] = int_(X,t)[T*.Qs]
% Parameters
    nP = 100 ; % source path discretization
    dh_dx = 1/4 ; % mesh density variation
    hmax = norm([L T])/30 ;
    interp = 'time' ; % 'time' or 'distance'
    distGeo = 'point' ; % 'line' or 'point'
% Domain
    lvlst = pkg.geometry.levelset.Rectangle([0 0 ; L T]) ;
% Source path discretization
    switch interp
        case 'distance'
            l = cumsum([0 ; sqrt(sum(diff(Xs).^2,2))]) ;
        case 'time'
            l = Xs(:,1) ;
    end
    lp = linspace(l(1),l(end),nP)' ;
    P = interp1(l,Xs,lp*.99) ;
% Mesh density function
    h0 = min(sqrt(min(sum(diff(P).^2,2))),hmax) ;
    switch distGeo
        case 'line'
            dist = @(p)pkg.geometry.distance.toPolyline(p,P) ;
        case 'point'
            Pp = permute(P,[3 2 1]) ;
            dist = @(p)min(pkg.geometry.distance.toPoint(p,Pp),[],3) ;
    end
    fh = @(p)min(h0 + dh_dx*dist(p),hmax) ;
% Fixed points
    pfix = [P ; lvlst.discretizeContour(fh)] ;
    pfix = uniquetol(pfix,3*h0/4,'DataScale',1,'ByRows',true) ;
% Build the mesh
    clf
    axis equal
    mesh_2D = lvlst.mesh('h0',h0 ...
                        ,'fh',fh ...
                        ,'pfix',pfix ...
                        ,'qmin',0 ...
                        ,'t_dmax',0 ...
                        ,'bndCons',false ...
                        ,'showMesh',true ...
                        ) ;
% Display
    clf reset
    axis equal, axis tight
    plot(mesh_2D,'FaceColor','none') ;
    xlabel 'Time (sec)'
    ylabel 'Position (m)'
    %plot(P(:,1),P(:,2),'.k','markersize',10)
    
%% SOLVE ON THE UNSTRUCTURED MESH
% Problem matrices
    [E,w,ie] = mesh_2D.integration ; W = spdiags(w(:),0,numel(w),numel(w)) ;
    N = mesh_2D.interpMat(E,ie) ;
    G = mesh_2D.gradMat(E,ie) ; Gy = G{2} ; Gt = G{1} ;
% Assemble (DO NOT SYMMETRIZE)
    K = S*lmbda*Gy'*W*Gy + S*Cp*rho*N'*W*Gt ;
% Second-hand term
    dt = [ts(2)-ts(1) ; (ts(3:end)-ts(1:end-2))/2 ; ts(end)-ts(end-1)] ;
    f = Qs*mesh_2D.interpMat(Xs)'*dt ;
% Boundary conditions
    keepDOF = mesh_2D.Nodes(:,1)>eps ;
    K = K(keepDOF,keepDOF) ;
    f = f(keepDOF) ;
% Solve
    temp_st = zeros(mesh_2D.nNodes,1) ;
    temp_st(keepDOF) = K\f ;
% Display
    clf reset
    pl = plot(mesh_2D) ;
    pl.FaceColor = 'interp' ;
    pl.Faces.FaceVertexCData = temp_st ;
    pl.VisibleEdges = 'none' ;
    pl.NodeCoordinates = mesh_2D.Nodes.*[1 1] ;
    %title 'Space-Time, Unstructured mesh'
    xlabel 'Time (sec)'
    ylabel 'Position (m)'
    colorbar
    axis equal,axis tight
    colormap(jet(16))
    caxis(Tlim)


    
    
    
    
    