%% TEST MESH OPERATIONS

%% CREATE A WIREFRAME MESH
    % Node coordinates: a sin curve
        x = linspace(0,1,10) ; x = [x(:) 0.1*cos(2*pi*x(:)) x(:)*0] ; % Sine curve
        %x = linspace(0,1*pi,10) ; x = 0.5*[cos(x(:)) sin(x(:)) x(:)*0] ; % Arc
    % Create the mesh automatically (wireframe mesh)
        mesh = pkg.geometry.mesh.Mesh(x) ;
    % Display
        cla ; axis equal ; pl = plot(mesh) ;
        xlabel x ; ylabel y ; zlabel z ;
        set(gca,'clipping','off') ;
        set(gca,'Interactions',[zoomInteraction rotateInteraction])

%% Offset the mesh with filling
    DIST = 0.05 ; fill = true ;
    mesh.offset(DIST,fill) ;
    pl.update ;

%% Extrude the mesh
    VEC = [0 0 1]*1 ; DIV = 50 ;
    mesh.extrude(VEC,DIV) ;
    pl.update ;
    
%% Sweep the mesh along a curve
    crv = linspace(0,1,51) ;
    crv = [crv(:).*0 0.1*cos(2*pi*crv(:)) crv(:)] ;
    mesh.sweep(crv,'angle') ;
    pl.update ; 
    
%% Mesh revolution
    AX = [1 0 0] ; PT = [0 -1 0] ; ANG = [0 pi] ; DIV = 20 ;
    mesh.revolution(AX,PT,ANG,DIV) ;
    pl.update ;
    
%% Mesh Slicing
    %lvlst = @(x)x(:,3) ; % horizontal plane at z=0
    %lvlst = @(x)x(:,3)-50/100 ; % horizontal plane at z=0.5
    lvlst = @(x)sum((x-[0 -10/10 1]).^2,2)-(9/10)^2 ; % sphere
    pl.Faces.FaceColor = 'flat' ; 
    pl.VisibleFaces = 'all' ; 
    pl.Faces.FaceVertexCData = mesh.featureLevelSetSign(lvlst,mesh.Faces) ;
    caxis([-1 1]) ;
    %%
    tic
    slice = mesh.cut(lvlst).OUT ;
    %slice = mesh.slice(lvlst) ;
    toc
    tag = 'SliceMesh' ; delete(findobj(gca,'tag',tag)) ;
    ppp = plot(slice,'Tag',tag)%,'HighlightEndNodes',true) ;
    
    
    
    
    
    