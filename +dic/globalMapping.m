function [N,INSIDE] = globalMapping(img,mesh,edgeMargin,tol) 
% Global mapping of a mesh in a function
% inputs:
%   - img is a 2D image
%   - mesh is a pkg.geometry.mesh.Mesh
%   - edgeMargin is used to map pixels outside the mesh (localized near edges)
%   - tol is the grid indices localization tolerance
% outputs:
%   - N is the SPARSE mapping matrix [nPixels nNodes]: shape functions evaluated
%       at valid pixel coordinates
%   - INSIDE is the SPARSE logical matrix [nPixels nElems] that gives the element associated to each pixel

% Default Parameters
    if nargin<3 ; edgeMargin = 0 ; end
    if nargin<4 ; tol = 1e-6 ; end 
    
% Infos (backup because the mesh may change)
    [nI,nJ] = size(img) ;
    nNodes = mesh.nNodes ;
    nElems = mesh.nElems ;
    
% edgeMargin: add quads as boundary edge outer offset
if edgeMargin
% Node normals
    [~,bndE,Nn,bndN] = mesh.boundaryNormals ;
% Keep only valid normals
    valid = ~any(isnan(Nn),2) ;
    Nn = Nn(valid,:) ;
    bndN(bndN) = valid ;
% Extrapolate normals
    newNodes = mesh.Nodes(bndN,:) + edgeMargin.*Nn(:,1:mesh.nCoord) ;
% Get edge indices
    newNodIdx = zeros(nNodes,1) ;
    newNodIdx(bndN) = (1:sum(bndN)) + nNodes ;
    edgNodIdx = mesh.Edges.subpart(bndE).NodeIdx ;
    offsetEdgIdx = reshape(newNodIdx(edgNodIdx),size(edgNodIdx)) ;
% New quads
    quadIdx = [flip(edgNodIdx,2) offsetEdgIdx] ;
% New Elems
    quadElems = pkg.geometry.mesh.elements.ElementTable(...
                            'Types',pkg.geometry.mesh.elements.base.Quadrangle ...
                            ,'Indices',padarray(quadIdx,[0 1],1,'pre') ...
                        ) ;
% Transfer matrix for new nodes
% Find CURRENT elements attached to CURRENT boundary nodes
    e2n = mesh.elem2node ;
    e2n = e2n(bndN,:) ;
    [nn,ee] = find(e2n) ;
% Localize newNodes with extrapolation
    [E,ie] = mesh.localize(newNodes ...
                            ,mesh.Elems ... localize in elements
                            ,true ... extrapolate
                            ,mesh.Nodes ... in the reference config
                            ,tol ... with reduced tolerance
                            ,nn,ee ... 
                            ) ;  
% Interpolation matrix T [nNodes+nNewNodes nNodes] (see below)
    T = mesh.interpMat(E,ie) ; % [nLocalizations nNodes]
    l2n = sparse(nn,(1:numel(nn))',1,size(newNodes,1),numel(nn)) ; % [nNewNodes nLocalizations]
    T = (l2n*T)./sum(l2n,2) ; % take the mean over localizations
    T = [sparse(1:nNodes,1:nNodes,1) ; T] ; % add old nodes
% Keep some element info before modification
    ele2edg = mesh.elem2edge ;
    [~,bndElem] = find(ele2edg(bndE,:)) ;
% Take a LOCAL copy of the mesh and modify it
    mesh = copy(mesh) ;
    mesh.Nodes = [mesh.Nodes ; newNodes] ;
    mesh.Elems = [mesh.Elems quadElems] ;
end

% Pixel localization
% ELEMENT bounding boxes
    xe = mesh.Elems.dataAtIndices(mesh.Nodes) ; % [nElems nMaxNodesByElem nCoord]
    bbox = [min(xe,[],2) max(xe,[],2)] ; % [nElems 2 nCoord]
    bbox = bbox + [-1 1]*tol ; % add the tolerance
% Which pixel is in which element bounding box ?
    [ji,ee] = pkg.data.domainIndices(bbox) ;
% Cull pixels out of image
    valid = all(ji>0 & ji<=[nJ nI],2) ;
    ji = ji(valid,:) ;
    ee = ee(valid) ;
% Localize in local element coordinates
    [E,ie] = mesh.localize(ji ...
                            ,mesh.Elems ... localize in elements
                            ,false ... do not extrapolate
                            ,mesh.Nodes ... in the reference config
                            ,tol ... with reduced tolerance
                            ,(1:numel(ee))' ... for all points
                            ,ee ... with respect to previously found elements
                            ) ;  
                        
% Cull invalid pixels
    valid = all(~isnan(E),2) ;
    E = E(valid,:) ;
    ie = ie(valid) ;
    ji = ji(valid,:) ;
    
% Linear indices
    idx = sub2ind([nI nJ],ji(:,2),ji(:,1)) ;  
    
% Sort linear indices
    [idx,is] = sort(idx,'ascend') ;
    E = E(is,:) ;
    ie = ie(is) ;
    
% Shape functions matrix
    Np = mesh.interpMat(E,ie) ;
    [pp,nn,vv] = find(Np) ;
    N = sparse(idx(pp),nn,vv,nI*nJ,mesh.nNodes) ;
    
% IF EDGE MARGINS HAVE BEEN ADDED
if edgeMargin 
% use transfer matrix for extrapolated nodes
    N = N*T ;
% Change the element number of pixel localized in boundary quad elements
    outIE = ie>nElems ;
    ie(outIE) = bndElem(ie(outIE)-nElems) ;
end
    
% INSIDE
    INSIDE = sparse(idx,ie,1,nI*nJ,nElems) ;
    
% Mean over localizations
    nLoc = sum(INSIDE,2) ;
    valid = find(nLoc>0) ;
    m = sparse(valid,valid,1./nLoc(valid),nI*nJ,nI*nJ) ;
    N = m*N ;
    
end