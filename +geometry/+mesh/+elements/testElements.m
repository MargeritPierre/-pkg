%% TEST CUSTOM ELEMENT

%% 1) CREATE THE ELEMENT AND PLOT IT
% Check node/edge/face numbering
    elmt = ... pkg.geometry.mesh.elements.base.Bar ...
           ... pkg.geometry.mesh.elements.base.Triangle ...
           ... pkg.geometry.mesh.elements.base.Quadrangle ...
           ... pkg.geometry.mesh.elements.base.Tetrahedron ...
           ... pkg.geometry.mesh.elements.base.Hexahedron ...
           ... pkg.geometry.mesh.elements.base.Prism ...
            pkg.geometry.mesh.elements.base.Pyramid ...
          ;

    clf
    axis equal
    plot(elmt) ;

%% 2) TEST IF POINTS ARE INSIDE THE REFERENCE ELEMENT 
% Check face/edge orientation

    nPts = 1000000 ; % number of test points
    extBbox = 0.2 ; % extend bbox
    tol = 1e-2 ; % to check if 'on' is working

    bbox = elmt.localCoordinatesDomain ;
    bbox = bbox + extBbox*range(bbox,1).*[-1; 1] ;
    
    E = (rand(nPts,elmt.nDims).*range(bbox,1)) + bbox(1,:) ;

    tic ; [in,on] = elmt.isInside(E,tol) ; toc
    
    disp(['Approximate Length/Area/Volume measurement: ' num2str(prod(range(bbox,1))*sum(in)/nPts)])
    
    clf
    axis equal
    plot(elmt) ;
    E = [E zeros(size(E,1),3-elmt.nDims)] ;
    plot3(E(in&~on,1),E(in&~on,2),E(in&~on,3),'.r','markersize',5) % STRICTLY INSIDE
    plot3(E(on,1),E(on,2),E(on,3),'.b','markersize',5) % ON BOUNDARY
    %plot3(E(~in,1),E(~in,2),E(~in,3),'.g','markersize',5) % OUTSIDE
    


%% 3) PLOT THE ELEMENT SHAPE FUNCTIONS AND LOW-ORDER DERIVATIVES
% Check the validity of shape functions and derivatives

    clf
    plot(elmt,'shapefunctions') ;
    
    
    
    
    
    