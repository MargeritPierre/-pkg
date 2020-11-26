% GENERATE A SERIES OF IMAGES CORRESPONDING TO A TEST

%% SPECIMEN GEOMETRY DEFINITION (LEVELSET)

%% DOG-BONE SPECIMEN
% Dimensions
    totalLength = 30 ;
    totalWidth = 20 ;
    usedWidth = 20 ;
    usedLength = 10 ;
    holeRadius = 2/10*usedWidth/2 ;
    
% Dependent dims
    R = (totalWidth-usedWidth)/2 ;
    L = (totalLength-usedLength)/2 ;
    bbox = [0 0 ; totalLength totalWidth] ;
    remRect = [L -R ; totalLength-L R] ;

% Define the shape
    lvlst = pkg.geometry.levelset.Rectangle(bbox) ;
    if usedWidth<totalWidth
        lvlst = lvlst - pkg.geometry.levelset.Rectangle(remRect) ;
        lvlst = lvlst - pkg.geometry.levelset.Rectangle(bbox(2,:)-remRect) ;
    end
    if R>0
        lvlst = lvlst - pkg.geometry.levelset.Circle([L 0],R) ;
        lvlst = lvlst - pkg.geometry.levelset.Circle([L totalWidth],R) ;
        lvlst = lvlst - pkg.geometry.levelset.Circle([totalLength-L 0],R) ;
        lvlst = lvlst - pkg.geometry.levelset.Circle([totalLength-L totalWidth],R) ;
    end
    if holeRadius>0
        lvlst = lvlst - pkg.geometry.levelset.Circle([totalLength totalWidth]/2,holeRadius) ;
    end
    
clf reset
axis equal, box on
h = plot(lvlst) ;


%% MESH THE SPECIMEN
delete(findobj(gcf,'type','hggroup')) ;

% Mesh edge length
    h0 = norm(range(lvlst.BoundingBox,1))/200 ;
% Apply distmesh
    mesh = lvlst.mesh('h0',h0) ;
% Display
    pl = plot(mesh) ;
    
    
%% GENERATE THE SPECIMEN MASK
    
% Parameters
    pixelsByMM = 40 ;
    margin = 20*[0 1] ; % IN PIXELS !
    
% Image size
    imgSz = pixelsByMM*flip(range(lvlst.BoundingBox,1)) + 2*flip(margin) ;
    imgBox = lvlst.BoundingBox + [-1;1].*margin/pixelsByMM ;
    
% Pixel coordinates
    px = linspace(imgBox(1,1),imgBox(2,1),imgSz(2)) ;
    py = linspace(imgBox(1,2),imgBox(2,2),imgSz(1)) ;
    [px,py] = meshgrid(px,py) ;
    P = [px(:) py(:)] ;

% MASK
    MASK = reshape(lvlst.inside(P),imgSz) ;
    
% Display
    clf
    axis ij
    axis equal
    axis tight
    image(repmat(MASK,[1 1 3])) ;
    
    
%% CREATE THE REFERENCE IMAGE
    
% Speckle image
    density = 7/10 ;
    medFiltSize = 3 ;
    gaussFiltSize = 7 ;
    edgeWidth = 2 ;
    backgroundLvl = 0.5 ;
    refImage = pkg.dic.generate.speckleImage(MASK,density,medFiltSize,gaussFiltSize,edgeWidth,backgroundLvl) ;
    
% Display
    clf
    axis ij
    axis equal
    axis tight
    image(repmat(refImage,[1 1 3])) ;
    
    
%% INTERPOLATION MATRIX

% Parameters
    edgeMargin = gaussFiltSize ;
    
% Convert Nodes in pixel coordinates
    interpMesh = copy(mesh) ;
    interpMesh.Nodes = (mesh.Nodes-imgBox(1,:))*pixelsByMM + 1  ;

% Localize pixels in the mesh
    profile on
    tic
    [N,INSIDE] = pkg.dic.globalMapping(MASK,interpMesh,edgeMargin) ;
    toc
    profile off
    
% Display
    clf
    axis ij
    axis equal
    axis tight
    if 0 % plot the sum Interpolation matrix (should be 1 in the domain)
        data = sum(N,2) ;
    elseif 1 % plot an interpolated function
        data =  N*interpMesh.Nodes(:,1) ...
               ... N*interpMesh.Nodes(:,2) ...
               ... N*u(:,1) ...
               ... N*u(:,2) ...
                ;
    else % plot the INSIDE mask
        data = sum(INSIDE,2) ;
    end
    I = reshape(full(data),size(MASK)) ;
    imagesc(I) ;
    pl = plot(interpMesh,'FaceColor','none') ;
    colorbar
    
    
%% GENERATE A DISPLACEMENT FIELD USING FEM
% Parameters
    nu = .3 ; % Poisson coeff.
    % Nodes clamped on the left
    clampedNodes = abs(mesh.Nodes(:,1)-min(mesh.Nodes(:,1)))<1e-6 ;
    % Nodes loaded on the right
    loadedNodes = abs(mesh.Nodes(:,1)-max(mesh.Nodes(:,1)))<1e-6 ;
    % Force the right edge to remain straight (uses e reference point)
    linkedNodes = loadedNodes ;
    uNorm = 0.005 ;
    
% Gauss integration objects
    [E,ie] = mesh.Elems.getListOf('GaussIntegrationPoints') ;
    nGP = numel(ie) ;
    W = mesh.Elems.getListOf('GaussIntegrationWeights') ;
    W = spdiags(W.*mesh.detJacobian(E,ie),0,nGP,nGP) ;
    
% Gradient matrix, df_dxi = G{i}*fn
    G = mesh.gradMat(E,ie) ;
% Strain matrix, eps = B*u
    O = sparse(nGP,mesh.nNodes) ;
    B = [G{1} O ; O G{2} ; G{2} G{1}] ;
% Stiffness Matrix
    C = [1 nu 0 ; nu 1 0 ; 0 0 (1-nu)/2] ;
    C = kron(sparse(C),W) ;
    K = B'*C*B ;
    
% Boundary load
    f = kron([1 ; 0],sparse(loadedNodes)) ; 
    
% Boundary Conditions & constraints using a transfer matrix
% T = [2*nNodes nDOF] at the end
    T = speye(2*mesh.nNodes) ; % [2*nNodes 2*nNodes]
    keepDOF = true(2*mesh.nNodes,1) ;
% Clamp
    T([clampedNodes ; clampedNodes],:) = 0 ;
    keepDOF([clampedNodes ; clampedNodes]) = false ;
% Reference points
    T([linkedNodes ; linkedNodes],:) = 0 ;
    T([linkedNodes ; false(mesh.nNodes,1)],end+1) = 1 ;
    T([false(mesh.nNodes,1) ; linkedNodes],end+1) = 1 ;
    keepDOF([linkedNodes ; linkedNodes]) = false ;
    keepDOF = [keepDOF ; true ; true] ;
% Cull unused DOF
    T = T(:,keepDOF) ;
% Apply BC
    K = T'*K*T ;
    f = T'*f ;

% Solve
    K = 0.5*(K + K') ; % force symmetry
    u = sparse(mesh.nNodes,mesh.nCoord) ;
    u(:) = T*(K\f) ;
    u = u.*uNorm./max(abs(u(:))).*norm(range(mesh.Nodes,1)) ;
    
% Strains
    EPS = reshape(B*u(:),nGP,3) ;
    EPS = mesh.interpMat(E,ie)\EPS ;
    
% Display
    clf
    axis equal
    dispMesh = copy(mesh) ; dispMesh.Nodes = dispMesh.Nodes+u ;
    pl = plot(dispMesh) ;    
    %pl.VisibleEdges = 'none' ;
    pl.Faces.FaceColor = 'interp' ;
    pl.Faces.FaceVertexCData = full(u(:,2))*pixelsByMM ; full(EPS(:,2)) ;
    %pl.Faces.FaceVertexCData = full(-EPS(:,2)./EPS(:,1)) ; caxis(nu+[-1 1]*0.01)
    colormap(jet(18))
    colorbar
    
    
%% DEFORMED IMAGE(S): FORWARD PROPAGATION (!!!! DISTORTIONS !!!)
% Parameters
    U = full(u).*reshape(linspace(0,1,20),1,1,[]) ; % [nNodes nCoord nIMG]
    
% Init images
    nIMG = size(U,3) ;
    IMG = zeros([imgSz nIMG],class(refImage)) ;
    
% Interpolation mesh
    interpMesh = copy(mesh) ;
    interpMesh.Nodes = (mesh.Nodes-imgBox(1,:))*pixelsByMM + 1  ;
    dU = cat(3,U(:,:,1),diff(U,1,3))*pixelsByMM ;
    
% Init Figure
    clf
    axis ij
    axis equal
    axis tight
    im = image(repmat(IMG(:,:,1),[1 1 3])) ;
    im.Interpolation = 'bilinear' ;
    pl = plot(interpMesh,'FaceColor','none') ;
    
% Genrate images
    for iii = 1:nIMG
        uu = N*U(:,:,iii) ;
        defImg = interp2(px,py,refImage,px(:)-uu(:,1),py(:)-uu(:,2),'linear',backgroundLvl) ;
        IMG(:,:,iii) = reshape(defImg,imgSz) ;
        im.CData = repmat(IMG(:,:,iii),[1 1 3]) ;
        interpMesh.Nodes = interpMesh.Nodes + dU(:,:,iii) ;
        pl.update ;
        drawnow ;
    end
    
    
%% DEFORMED IMAGE(S): EXACT INVERSE PROPAGATION
% Parameters
% Displacement U in pixels [nNodes nCoord nIMG]
    U = pixelsByMM*full(u).*reshape(linspace(0,1,2),1,1,[]) ;
    
% Init images
    nIMG = size(U,3) ;
    IMG = zeros([imgSz nIMG],class(refImage)) ;
    
% Interpolation mesh
    interpMesh = copy(mesh) ;
    interpMesh.Nodes = (mesh.Nodes-imgBox(1,:))*pixelsByMM + 1  ;
    dU = cat(3,U(:,:,1),diff(U,1,3)) ;
    
% Init Figure
    clf reset
    axis ij
    axis equal
    axis tight
    im = image(repmat(IMG(:,:,1),[1 1 3])) ;
    im.Interpolation = 'bilinear' ;
    pl = plot(interpMesh,'FaceColor','none') ;
    
% Genrate images
    for iii = 1:nIMG
        interpMesh.Nodes = interpMesh.Nodes + dU(:,:,iii) ;
        [Ni,INi] = pkg.dic.globalMapping(MASK,interpMesh,edgeMargin) ;
        X = [px(:) py(:)] - Ni*U(:,:,iii)/pixelsByMM ;
        defImg = interp2(px,py,refImage,X(:,1),X(:,2),'cubic',backgroundLvl) ;
        defImg(sum(INi,2)==0) = backgroundLvl ;
        IMG(:,:,iii) = reshape(defImg,imgSz) ;
        im.CData = ... reshape(full(sum(Ni,2)),imgSz) ...
                    repmat(IMG(:,:,iii),[1 1 3]) ...
                   ;
        pl.update ;
        drawnow ;
    end
    
%% IMAGE SLIDER
    
    clf reset
    axis ij
    axis equal
    axis tight
    im = image(repmat(IMG(:,:,1),[1 1 3])) ;
    im.Interpolation = 'bilinear' ;
    
    slider = uicontrol(gcf,'Style','slider','Min',1,'Max',nIMG,'Value',1) ;
    slider.Units = 'normalized' ;
    slider.Position = [0 0 1 0.04] + [1 1 -2 -2]*0.001 ;
    addlistener(slider,'Value','PostSet',@(src,evt)set(im,'CData',repmat(IMG(:,:,round(slider.Value)),[1 1 3])))

    
%% SAVE IMAGES

    [file,path] = uiputfile('img_%03i.png','Save the images') ;
    if path==0 ; return ; end
    
    for iii = 1:nIMG
        filename = [path sprintf(file,iii)] ;
        imwrite(IMG(:,:,iii),filename) ;
    end
    

    

    
    

    
    













