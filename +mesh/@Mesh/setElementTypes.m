function mesh = setElementTypes(this,types)
% Change the types of elements
% Return a copy of the mesh if an output is queried (the original mesh will
% remain unchanged)
% <TODO> how MeshFunctions can react to this to stay valid ?
if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end

% Get the current element types
    oldElmtTypes = mesh.Elems.Types ;
    nElmtTypes = numel(mesh.Elems.Types) ;
    
% Format the current type list
    if numel(types)==1
        types = repmat(types,[nElmtTypes 1]) ;
    elseif numel(types)~=numel(oldElmtTypes)
        error('he number of provided new element types must match the current one.')
    end
    newElmtTypes = types ;
    
% Interpolation matrices for nodal values
    iMat = {} ;
    for ee = 1:nElmtTypes
        iMat{ee} = oldElmtTypes(ee).evalAt(newElmtTypes(ee).NodeLocalCoordinates) ;
    end
    
% New list of nodes coordinates and indices
    Xe = mesh.Elems.dataAtIndices(mesh.X.Values) ; % [nElems nMaxNodesByElem nCoord] ;
    X = [] ;
    NodeIdx = zeros(mesh.nElems,max([newElmtTypes.nNodes]),'uint32') ;
    for ee = 1:nElmtTypes
        % Elements with the given type
            ie = find(mesh.Elems.TypeIdx==ee) ; 
            nElmts = numel(ie) ;
        % New/Old element infos
            nOldNodes = oldElmtTypes(ee).nNodes ;
            nNewNodes = newElmtTypes(ee).nNodes ;
        % Corresponding coordinates at element nodes
            Xee = Xe(ie,1:nOldNodes,:) ; % [nEi nOldNodes nCoord] ;
            Xee = permute(Xee,[2 1 3]) ; % [nOldNodes nEi nCoord] ;
            nXee = reshape(iMat{ee}*reshape(Xee,nOldNodes,[]),nNewNodes,nElmts,[]) ; % [nNewNodes nEi nCoord] ;
            nXee = permute(nXee,[2 1 3]) ; % [nEi nNewNodes nCoord] ;
        % New Node Indices
            NodeIdx(ie,1:nNewNodes) = size(X,1) + reshape(1:numel(nXee(:,:,1)),size(nXee(:,:,1))) ;
            X = [X ; reshape(nXee,[],mesh.nCoord)] ;
    end

% Cull duplicated nodes
    [X,~,nn] = uniquetol(X,mesh.defaultTolerance,'ByRows',true,'DataScale',1) ;
    NodeIdx = reshape(nn(NodeIdx),mesh.nElems,max([newElmtTypes.nNodes])) ;
    
% Assign to the element
    mesh.X.Values = X ;
    mesh.Elems = pkg.mesh.elements.ElementTable('Types',newElmtTypes,'Indices',[mesh.Elems.TypeIdx NodeIdx]) ;
    

end

