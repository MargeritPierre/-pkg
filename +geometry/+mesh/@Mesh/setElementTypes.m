function this = setElementTypes(this,types)
% Change the types of elements
% Return a copy of the mesh if an output is queried (the original mesh will
% remain unchanged)
% <TODO> how MeshFunctions can react to this to stay valid ?

% Get the current element types
    oldElmtTypes = this.Elems.Types ;
    nElmtTypes = numel(this.Elems.Types) ;
    
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
    Xe = this.Elems.dataAtIndices(this.Nodes) ; % [nElems nMaxNodesByElem nCoord] ;
    X = [] ;
    NodeIdx = zeros(this.nElems,max([newElmtTypes.nNodes]),'uint32') ;
    for ee = 1:nElmtTypes
        % Elements with the given type
            ie = find(this.Elems.TypeIdx==ee) ; 
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
            X = [X ; reshape(nXee,[],this.nCoord)] ;
    end
    
% Create the new mesh
    mesh = copy(this) ;
    mesh.Nodes = X ;
    mesh.Elems = pkg.geometry.mesh.elements.ElementTable('Types',newElmtTypes,'Indices',[this.Elems.TypeIdx NodeIdx]) ;

% Cull duplicated nodes without changing the element-element connectivity
    % Get nodes that are duplicated in the new mesh
        [X,un,nn] = uniquetol(mesh.Nodes,mesh.defaultTolerance,'ByRows',true,'DataScale',1) ;
    % welding Matrix
        W = sparse(nn,1:numel(nn),1,numel(un),numel(nn)) ;
    % Find the element-element connectivity changes
        old_el2no = this.elem2node ;
        new_el2no = W*mesh.elem2node ;
        old_el2el = logical(old_el2no'*old_el2no)-speye(this.nElems) ;
        new_el2el = logical(new_el2no'*new_el2no)-speye(mesh.nElems) ;
        [el1,el2] = find(logical(old_el2el-new_el2el)) ;
    % Get associated elements couples
        el = [el1(:) el2(:)] ;
        el = el(el(:,1)~=el(:,2),:) ;
        el = unique(sort(el,2),'rows') ;
    
% % Restablish special connectivities
%     % Find the welded nodes
%         [welded,~] = find(el2no(:,el(:,1)) & el2no(:,el(:,2))) ;
%     % Split the welded nodes
%         mesh.splitNodes(unique(welded)) ;
        
    X = W.*(1./sum(W,2))*X ;
    NodeIdx = reshape(nn(NodeIdx),this.nElems,[]) ;

% % Cull duplicated nodes
% % /!\ this is welding mesh features that were splitted !
%     [X,~,nn] = uniquetol(X,this.defaultTolerance,'ByRows',true,'DataScale',1) ;
%     NodeIdx = reshape(nn(NodeIdx),this.nElems,[]) ;
    
% % Create the new mesh
%     mesh = copy(this) ;
%     mesh.Nodes = X ;
%     mesh.Elems = pkg.geometry.mesh.elements.ElementTable('Types',newElmtTypes,'Indices',[this.Elems.TypeIdx NodeIdx]) ;

% Replace this ?
    if nargout==0 
        this.Nodes = mesh.Nodes ;
        this.Elems = mesh.Elems ;
    end

end

