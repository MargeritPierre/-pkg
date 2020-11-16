%% GENERATE A SERIES OF IMAGES CORRESPONDING TO A TEST

%% SPECIMEN GEOMETRY DEFINITION (LEVELSET)

%% DOG-BONE SPECIMEN
% Dimensions
    totalLength = 40 ;
    totalWidth = 20 ;
    usedWidth = 10 ;
    usedLength = 20 ;
    
% Dependent dims
    R = (totalWidth-usedWidth)/2 ;
    L = (totalLength-usedLength)/2 ;
    bbox = [0 0 ; totalLength totalWidth] ;
    remRect = [L -R ; totalLength-L R] ;

% Define the shape
    lvlst = pkg.geometry.levelset.Rectangle(bbox) ;
    lvlst = lvlst - pkg.geometry.levelset.Rectangle(remRect) ;
    lvlst = lvlst - pkg.geometry.levelset.Rectangle(bbox(2,:)-remRect) ;
    lvlst = lvlst - pkg.geometry.levelset.Circle([L 0],R) ;
    lvlst = lvlst - pkg.geometry.levelset.Circle([L totalWidth],R) ;
    lvlst = lvlst - pkg.geometry.levelset.Circle([totalLength-L 0],R) ;
    lvlst = lvlst - pkg.geometry.levelset.Circle([totalLength-L totalWidth],R) ;

clf
axis equal, box on
h = plot(lvlst) ;


%% MESH THE SPECIMEN
    h0 = norm(range(lvlst.BoundingBox,1))/200 ;
    delete(findobj(gcf,'type','hggroup')) ;
    mesh = lvlst.mesh('h0',h0,'debug',true,'bndCons',false) ;
    plot(mesh) ;


