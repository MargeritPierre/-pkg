%% TEST MESH PACKAGE


%% THE RANDOM MESH (code timing)

    clc,clearvars
    import pkg.mesh.*
    import pkg.mesh.elements.*
    
    nNodes = 10000000 ; % Number of nodes
    nElems = 10000 ; % Number of random elements
    elemType = LagrangeElement('quad',1) ; % Element type (only on type supported here)
    
    nodeIndices = randi([1 nNodes],[nElems elemType.nNodes]) ; % Indices of nodes connected to each element
    indices = [ones(size(nodeIndices,1),1) nodeIndices] ; % Including the element type


    tic ;
    %profile on
    
    % Create the element table
    elemTable = ElementTable('Types',elemType,'Indices',indices) ;
    
    % Create the mesh
    mesh = Mesh('Elems',elemTable) ;

    %profile viewer
    toc
    
    
%% THE MIXED-ELEMENT MESH (TRIS & QUADS)
    clc,clearvars
    import pkg.mesh.*
    import pkg.mesh.elements.*
    
    elemTypes = [LagrangeElement('tri',1) LagrangeElement('quad',1)] ;
    nodeIdx = [1 2 4 3 ; 3 7 5 NaN ; 4 3 6 5] ;
    typeIdx = [2 ; 1 ; 2] ;
    
    if 0 % Give the element type
        mesh = Mesh('Elems',ElementTable('Types',elemTypes,'Indices',[typeIdx nodeIdx])) ;
    else % Element type automatically assigned (verify with mesh.Elems.TypeIdx==typeIdx)
        mesh = Mesh('Elems',ElementTable('Types',elemTypes,'NodeIdx',nodeIdx)) ;
    end
    
    
    
    
    

