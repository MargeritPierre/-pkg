function mesh = setElementOrder(this,order)
% Change the order of LagrangeElements in the mesh
% Order is a positive or null SCALAR integer  
% Return a copy of the mesh if an output is queried (the original mesh will
% remain unchanged)
% <TODO> how MeshFunctions can react to this to stay valid ?
if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end

% Get the current element types
    oldElmtTypes = mesh.Elems.Types ;
    nElmtTypes = numel(mesh.Elems.Types) ;

% Depending on the mesh element class
    switch class(oldElmtTypes)
        case 'pkg.mesh.elements.LagrangeElement'
            % Check if the order shift is possible
                if range(oldElmtTypes.Order)>0 ; error('Order change on mixed-order Lagrange element meshes is not supported') ; end
            % Create the element types of replacement
                newElmtTypes = pkg.mesh.elements.AbstractElement.empty ;
                for ee = 1:nElmtTypes
                    newElmtTypes(end+1) = pkg.mesh.elements.LagrangeElement(oldElmtTypes.Geometry,order) ;
                end
        otherwise
            error('The mesh element types are not supported for order change') ;
    end
    
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
        ie = find(mesh.Elems.TypeIdx==ee) ; 
        nOldNodes = oldElmtTypes(ee).nNodes ;
        nNewNodes = newElmtTypes(ee).nNodes ;
        nElmts = numel(ie) ;
        Xee = Xe(ie,1:nOldNodes,:) ; % [nEi nOldNodes nCoord] ;
        Xee = permute(Xee,[2 1 3]) ; % [nOldNodes nEi nCoord] ;
        nXee = reshape(iMat{ee}*reshape(Xee,nOldNodes,[]),nNewNodes,nElmts,[]) ; % [nNewNodes nEi nCoord] ;
        nXee = permute(nXee,[2 1 3]) ; % [nEi nNewNodes nCoord] ;
        NodeIdx(ie,1:nNewNodes) = size(X,1) + reshape(1:numel(nXee(:,:,1)),size(nXee(:,:,1))) ;
        X = [X ; reshape(nXee,[],mesh.nCoord)] ;
    end
    [X,~,nn] = uniquetol(X,mesh.defaultTolerance,'ByRows',true,'DataScale',1) ;
    NodeIdx = reshape(nn(NodeIdx),mesh.nElems,max([newElmtTypes.nNodes])) ;
    
% Assign to the element
    mesh.X.Values = X ;
    mesh.Elems = pkg.mesh.elements.ElementTable('Types',newElmtTypes,'Indices',[mesh.Elems.TypeIdx NodeIdx]) ;
    

end

