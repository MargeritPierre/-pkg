function [ival,iint] = inDomain(values,domains,method)
%INDOMAIN return indices of nd-values contained in nd-domains
% values: [nVal nCoord]
% domains [nInt 2 nCoord] (minValue_d maxValue_d)
%
% Alternative implementation to the following "brute-force" test:
%   dom = permute(domains,[2 1 3]) ; % [2 nInt nCoord]
%   pt = reshape(values,[nVal 1 nCoord]) ; % [nVal 1 nCoord]
%   INSIDE = all(dom(1,:,:)<=pt,3) & all(dom(2,:,:)>=pt,3) ;
%	[ip,ie] = find(INSIDE) ;
%
% It is memory-efficient and the speedup depends on the sparsity of INSIDE
%   high sparsity: up to two order of magnitude faster!
%   low sparsity: hem.. not so efficient!
%

% Infos
nVal = size(values,1) ;
nInt  = size(domains,1) ;
nCoord = min(size(values,2),size(domains,3)) ;

% Format
if size(values,2)>nCoord ; values = values(:,1:nCoord) ; end
if size(domains,3)>nCoord ; domains = domains(:,1:2,1:nCoord) ; end

% Choose the best method
if nargin<3 
    if nVal<=1 || nInt<=1
        method = 'brute' ;
    elseif nCoord<=1 % 1D domains
        if nVal*nInt<1e7
            method = 'brute' ;
        else % sorting is perfect for 1D problems
            method = 'sort' ; 
        end
    elseif nCoord<=1 || nVal<nInt
        method = 'sort' ; 
    else
        method = 'tree' ;
    end
end

% Process
switch method
    case 'sort'
        coordSorting ;
    case 'tree'
        kdTree ;
    case 'tree2'
        extTree ;
    otherwise
        bruteForce ;
end


% END OF MAIN FUNCTION
% ===================================================


%% NESTED FUNCTIONS 


% BRUTE FORCE ALGORITHM: TEST ALL vs. ALL
function bruteForce
    if nVal*nInt<=1e6 % brute fully vectorized
        dom = permute(domains,[2 1 3]) ; % [2 nInt nCoord]
        pt = reshape(values,[nVal 1 nCoord]) ; % [nVal 1 nCoord]
        INSIDE = all(dom(1,:,:)<=pt,3) & all(dom(2,:,:)>=pt,3) ;
        [ival,iint] = find(INSIDE) ;
    else
        dom = permute(domains,[1 3 2]) ; % [nInt nCoord 2]
        pt = values ; % [nVal nCoord]
        if nVal<nInt % test all domains vs. one value at once
            iint = cell(nVal,1) ;
            dmin = dom(:,:,1) ;
            dmax = dom(:,:,2) ;
            for vv = 1:nVal
                vvv = pt(vv,:) ;
                in = all(dmin<=vvv & dmax>=vvv,2) ;
                iint{vv} = find(in) ;
            end
            ni = cellfun(@numel,iint) ;
            ival = repelem((1:nVal)',ni,1) ;
            iint = cat(1,iint{:}) ;
        else % test all values vs. one domain at once
            ival = cell(nInt,1) ;
            for ii = 1:nInt
                ddd = dom(ii,:,:) ;
                in = all(ddd(:,:,1)<=pt & ddd(:,:,2)>=pt,2) ;
                ival{ii} = find(in) ;
            end
            nv = cellfun(@numel,ival) ;
            iint = repelem((1:nInt)',nv,1) ;
            ival = cat(1,ival{:}) ;
        end
    end
end


% FIND POINTS IN DOMAIN BY FIRST SORTING (MEMORY EFFICIENT)
% critical point: array accumulation (sparse)
function coordSorting

    % Put values and interval bounds together
    data = [values ; reshape(domains,[2*nInt nCoord])] ; % [nVal+2*nInt nCoord]
    nData = nVal+2*nInt ;
    % Keep indices corresponding to interval bounds
    ibounds = reshape(nVal+1:nData,[nInt 2]) ; % [nInt 2]

    % Sort data
    [~,isort] = sort(data,1) ; % [nData nCoord]
    % Retrieve indices corresponding to values
    isVal = isort<=nVal ; % [nData nCoord]

    % Number of values before each index
    np = cumsum(isVal(:)) ; % [nData*nCoord 1]
    % Sorted index of each interval bound
    ii = NaN(nInt,2,nCoord) ;
    for cc = 1:nCoord ; [~,ii(:,:,cc)] = ismember(ibounds,isort(:,cc)) ; end
    ii = ii + reshape((2*nInt+nVal)*(0:nCoord-1),[1 1 nCoord]) ; % [nInt 2 nCoord]
    % Number of values before the interval bound
    npi = np(ii) ; % [nInt 2 nCoord]
    npi = permute(npi,[1 3 2]) ; % [nInt nCoord 2]

    % Long vector of interval indices
    nIndicesInInterval = diff(npi,1,3) ; % [nInt nCoord]
    iint = repelem(repmat(1:nInt,[1 nCoord]),nIndicesInInterval(:)) ; % a lot of elements !

    % Long vector of value indices
    starts = npi(:,:,1)+1 ; ends = npi(:,:,2) ; % [nInt nCoord]
    lengths = ends(:)'-starts(:)'+1 ; cumLengths = cumsum(lengths) ; % [nInt*nCoord]
    vv = (0:cumLengths(end)-1) + repelem(starts(:)' - [0 cumLengths(1:end-1)],lengths) ; % lot of elements ! (critical assignment)

    % Retrieve the corresponding value indices
    ival = isort(isVal) ; % [nVal*nCoord]
    ival = ival(vv) ; % lot of elements ! (critical assignment)

    if nCoord>1 % multiple coordinates: accumulate repeated [val int] couples
        % Sort linear indices
        ii = sub2ind([nVal nInt],ival(:),iint(:)) ;
        ii = sort(ii) ;
        % Find places where indices change
        change = find([true ; ii(1:end-1)~=ii(2:end) ; true]) ;
        % Is there nCoord indices ?
        cc = diff(change)==nCoord ;
        % Keep only valid subindices
        [ival,iint] = ind2sub([nVal nInt],ii(change(cc))) ;
    end

end

% USING A KD-TREE FOR THE BOXES
function kdTree
    
    % Permute the domains
    dom = permute(domains,[1 3 2]) ; %[nInt nCoord 2) ;
    
    % Build the kd-tree
    div = 10 ; 
    order = max(floor(log(nInt)/log(div)) - 1,0) ;
    tree = pkg.graph.tree.KD(mean(dom,3),div,order) ;
    
    % Compute enclosing boxes corresp. to each tree node
    encBox = @(ii)cat(3,min(dom(ii,:,1),[],1),max(dom(ii,:,2),[],1)) ;
    encBoxes = cellfun(encBox,{tree.Nodes.Data},'UniformOutput',false) ;
    
    % Add the last layer with one node by domain
    en = find(tree.endNodes) ;
    dd = {tree.Nodes(en).Data} ;
    nd = cellfun(@numel,dd) ;
    newEdges = [reshape(repelem(en(:),nd),[],1) cat(2,dd{:})'+tree.nNodes] ;
    tree.Edges = [tree.Edges ; newEdges] ;
    tree.Nodes(end+nInt) = pkg.graph.Node ;
    
    % fill the tree with bounding boxes
    [tree.Nodes(1:end-nInt).Data] = deal(encBoxes{:}) ;
    cellDom = mat2cell(dom,ones(1,nInt)) ;
    [tree.Nodes(end-nInt+1:end).Data] = deal(cellDom{:}) ;
    
    % Dispatch the test in the tree
    allBoxes = cat(1,tree.Nodes.Data) ;
    nodeLayer = tree.nodeLayer ;
    parentNode = tree.parentNode ;
    IDX = cell(1,tree.nNodes) ;
    for ll = 1:max(nodeLayer)
    % nodes of the current layer
        nn = find(nodeLayer==ll) ;
    % corresp. parent node indices
        if ll==1
            pp = repmat({(1:nVal)'},[1 numel(nn)]) ;
        else
            pp = IDX(parentNode(nn)) ;
        end
    % flatten into a common list
        np = cellfun(@numel,pp) ; % number of indices in each parent node
        ppp = cat(1,pp{:}) ;
    % retrieve corresp.values and intervals
        vvv = values(ppp,:) ;
        dmin = repelem(allBoxes(nn,:,1),np,1) ;
        dmax = repelem(allBoxes(nn,:,2),np,1) ;
    % test if point in box
        inside = vvv>=dmin & vvv<=dmax ;
        inside = all(inside,2) ;
        if ~any(inside) ; break ; end
    % re-put in a cell array
        iii = reshape(repelem(1:numel(nn),np),[],1) ;
        ii = accumarray(iii(inside),ppp(inside),[numel(nn) 1],@(x){x(:)},{}) ;
    % record
        [IDX{nn}] = deal(ii{:}) ;
    end
    
    % Convert to indices
    pp = IDX(tree.endNodes) ;
    np = cellfun(@numel,pp) ;
    iint = repelem((1:nInt)',np) ;
    ival = cat(1,pp{:}) ;
    
end

% USING A KD-TREE FOR THE BOXES
function extTree
    
    % Add the index in the lists
    % This avoid critical indexing steps
    val = [(1:nVal)' values] ;
    
    % Build the kd-tree
    div = 10 ; % list division at each tree layer
    order = max(floor(log(nInt)/log(div)) - 1,0) ; % number of layers
    tree = pkg.graph.tree.KD(reshape(mean(domains,2),[nInt nCoord]),div,order) ;
    
    % Permute the data (memory access optimization)
    dom = permute(domains,[3 1 2]) ; % [nCoord nInt 2] ;
    val = permute(val,[2 1]) ; % [nCoord nVal]
    
    % Compute enclosing boxes corresp. to each tree node
    encBox = @(ii)cat(3,min(dom(:,ii,1),[],2),max(dom(:,ii,2),[],2)) ;
    encBoxes = cellfun(encBox,{tree.Nodes.Data},'UniformOutput',false) ;
    
    % Add the last layer with one node by domain
    en = find(tree.endNodes) ;
    dd = {tree.Nodes(en).Data} ;
    nd = cellfun(@numel,dd) ;
    newEdges = [reshape(repelem(en(:),nd),[],1) cat(2,dd{:})'+tree.nNodes] ;
    tree.Edges = [tree.Edges ; newEdges] ;
    
    % Dispatch the test in the tree
    allBoxes = [cat(2,encBoxes{:}) dom] ;
    nodeLayer = tree.nodeLayer ;
    parentNode = tree.parentNode ;
    IDX = cell(1,tree.nNodes) ;
    for ll = 1:max(nodeLayer)
    % nodes of the current layer
        nn = find(nodeLayer==ll) ;
    % corresp. parent node indices
        if ll==1
            pp = repmat({val},[1 numel(nn)]) ;
        else
            pp = IDX(parentNode(nn)) ;
        end
    % flatten into a common list
        [~,np] = cellfun(@size,pp) ; % number of values in each parent node
        ppp = cat(2,pp{:}) ;
    % retrieve corresp.values and intervals
        vvv = ppp(2:end,:) ;
        dmin = repelem(allBoxes(:,nn,1),1,np) ;
        dmax = repelem(allBoxes(:,nn,2),1,np) ;
    % test if point in box
        inside = all(vvv>=dmin & vvv<=dmax,1) ;
        if ~any(inside) ; break ; end
    % Count the number of insides in each node
        nInside = [0 cumsum(inside)] ;
        nInside = nInside(cumsum(np)+1) ;
        nInside = diff([0 nInside]) ;
    % re-put in a cell array
        ppp = ppp(:,inside) ;
    % record
        ends = cumsum(nInside) ;
        starts = [0 ends(1:end-1)] + 1 ;
        for nnn = find(nInside) ; IDX{nn(nnn)} = ppp(:,starts(nnn):ends(nnn)) ; end
    end
    
    % Convert to indices
    pp = IDX(end-nInt+1:end) ;
    [~,np] = cellfun(@size,pp) ;
    iint = repelem((1:nInt)',np) ;
    ival = cat(2,pp{:}) ;
    if~isempty(ival) ; ival = ival(1,:)' ; end
    
end

end
% END OF CODE
% ===================================================


function tests
%% TESTS ===================================================
clearvars
nC = 3 ;
nP = 1000000 ;
nB = 1000 ; Bsz = .1*(1/nB)^(1/nC) ;

if 0 % random boxes
    P = rand(nP,nC) ;
else % uniform grid
    nPd = ceil(nP^(1/nC)) ; 
    x = repmat({linspace(0,1,nPd)},[1 nC]) ; 
    [x{:}] = ndgrid(x{:}) ;
    P = reshape(cat(nC+1,x{:}),[nPd^nC nC]) ;
    nP = size(P,1) ;
end
if 0 % random boxes
    B = rand(nB,1,nC) + Bsz*[-1 1].*rand(nB,1,nC) ;
elseif 1 % box grid
    nBd = ceil(nB^(1/nC)) ; 
    x = repmat({linspace(0,1-1/nBd,nBd)},[1 nC]) ; 
    [x{:}] = ndgrid(x{:}) ;
    x = reshape(cat(nC+1,x{:}),[nBd^nC 1 nC]) ;
    B = [x x+1/nBd] ;
    nB = size(B,1) ;
end


if 0 && nB<1e5 && nP<1e5 && nC==2 % nC>1 && nC<4
    clf
    axis equal
    patch('vertices',P,'faces',(1:nP)','marker','.','MarkerEdgeColor','r')
    ind = permute(B,[1 3 2]) ;
    ind = [ind(:,:,1) ; [ind(:,1,2) ind(:,2,1)] ; ind(:,:,2) ; [ind(:,1,1) ind(:,2,2)]] ;
    patch('vertices',ind,'faces',(1:nB)'+(0:3)*nB,'FaceColor','none','EdgeColor','b')
end


clc
N = log10(nP*nB)
profile on
%tic ; [ival1,iint1] = pkg.data.inDomain(P,B,'brute') ; toc
%tic ; [ival2,iint2] = pkg.data.inDomain(P,B,'sort') ; toc
tic ; [ival3,iint3] = pkg.data.inDomain(P,B,'tree') ; toc
tic ; [ival4,iint4] = pkg.data.inDomain(P,B,'tree2') ; toc
profile off

%numel(iint2)
%nnz(sparse(ival1,iint1,true,nP,nB)-sparse(ival2,iint2,true,nP,nB))
%nnz(sparse(ival1,iint1,true,nP,nB)-sparse(ival3,iint3,true,nP,nB))
%nnz(sparse(ival1,iint1,true,nP,nB)-sparse(ival4,iint4,true,nP,nB))
%nnz(sparse(ival3,iint3,true,nP,nB)-sparse(ival4,iint4,true,nP,nB))



%% PLOT GROUPS
clf ;
axis equal
for gg = 1:numel(groups)
    ii = groups{gg} ;
    isValue = ii<=nVal ;
    ival = ii(isValue) ;
    med = plot(values(ival,1),values(ival,2),'.') ;
    iint = mod(ii(~isValue)-nVal-1,nInt)+1 ;
    iint = unique(iint) ;
    ind = permute(domains(iint,:,:),[1 3 2]) ;
    ind = [ind(:,:,1) ; [ind(:,1,2) ind(:,2,1)] ; ind(:,:,2) ; [ind(:,1,1) ind(:,2,2)]] ;
    patch('vertices',ind,'faces',(1:numel(iint))'+(0:3)*numel(iint),'FaceColor','none','EdgeColor',med.Color)
end

%% SPARSITY INDICATORS
minOvlp = reshape(sum(diff(B,1,2),1),[],nC)./range(reshape(B,[],nC),1) ;
sparsity = sum(prod(diff(B,1,2),3))/(prod(range(reshape(B,[],nC),1),2)*nB) ;





%% END OF TESTS
end





