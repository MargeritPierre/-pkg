function tree = KD(data,div,order,splitSpace)
%KD Build a kd tree from a point cloud
% data: [nData nCoord]
% div: number of partitions of the data at each recursion (median==2)
% order: tree.nLayers = order+1
% splitSpace: "index" or "coordinate"

if nargin==0 ; error('Not enough arguments') ; end

% Infos
data = data(:,:) ; % [nData nCoord]
[nData,nCoord] = size(data) ;

% Default division: 2 (median) 
if nargin<2 ; div = 2 ; end
% Automatic order: maximum order tree
if nargin<3 ; order = max( ceil( log(nData)/log(div) ) -  1 ,0) ; end
if nargin<4 ; splitSpace = 'coordinate' ; end

% Tree nodes
    % number of nodes by tree layer
    nNodesInLayer = div.^(0:order) ;
    % Total of nodes in the tree
    nNodesTotal = sum(nNodesInLayer) ;

% Tree edges
if nNodesTotal>1
    % previous node in the tree
    prevNode = [0 ; repelem(1:nNodesTotal-div^order,div)'] ;
    % edges (removing the root node ingoing edge)
    edges = [prevNode(2:end) (2:nNodesTotal)'] ;
else
    edges = [] ;
end

% Data splitting
    % Data mapping in [0 1[
    data = (data-min(data,[],1))./range(data,1)*(1-1/nData) ;
    % node indices: [nData nLayers]
    nodeIdx = [ones(nData,1) zeros(nData,order)] ;
    for oo = 1:order
        coord = mod(oo-1,nCoord)+1 ;
        ni = (nodeIdx(:,oo)-min(nodeIdx(:,oo),[],1) + data(:,coord))*div ;
        switch splitSpace
            case 'index'
                dd = nodeIdx(:,oo) + data(:,coord) ;
                [~,si] = sort(dd) ;
                ni = floor(linspace(0,1-1/nData,nData)*(div^oo)) + sum(nNodesInLayer(1:oo)) + 1 ;
                nodeIdx(si,oo+1) = ni ;
            case 'coordinate'
                nodeIdx(:,oo+1) = floor(ni) + sum(nNodesInLayer(1:oo)) + 1 ;
        end
    end
    nDataInNode = accumarray(nodeIdx(:),1,[nNodesTotal 1]) ;
    [~,dataIdx] = sort(nodeIdx,1) ;
    IDX = mat2cell(dataIdx(:),nDataInNode) ;

% Data splitting
    % Data mapping in [0 1[
    data = (data-min(data,[],1))./range(data,1)*(1-1/nData) ;
    % node indices: [nData nLayers]
    nodeIdx = linspace(0,1-1/nData,nData)'*nNodesInLayer(:)' ;
    nodeIdx = [ones(nData,1) zeros(nData,order)] ;
    for oo = 1:order
        coord = mod(oo-1,nCoord)+1 ;
        ni = (nodeIdx(:,oo)-min(nodeIdx(:,oo),[],1) + data(:,coord))*div ;
        switch splitSpace
            case 'index'
                dd = nodeIdx(:,oo) + data(:,coord) ;
                [~,si] = sort(dd) ;
                ni = floor(linspace(0,1-1/nData,nData)*(div^oo)) + sum(nNodesInLayer(1:oo)) + 1 ;
                nodeIdx(si,oo+1) = ni ;
            case 'coordinate'
                nodeIdx(:,oo+1) = floor(ni) + sum(nNodesInLayer(1:oo)) + 1 ;
        end
    end
    nDataInNode = accumarray(nodeIdx(:),1,[nNodesTotal 1]) ;
    [~,dataIdx] = sort(nodeIdx,1) ;
    IDX = mat2cell(dataIdx(:),nDataInNode) ;

% Build the tree
    tree = pkg.graph.tree.Tree ;
    % Edges
    tree.Edges = edges ;
    % Node data
    tree.Nodes = repelem(pkg.graph.Node,nNodesTotal) ;
    [tree.Nodes.Data] = deal(IDX{:}) ;



end

%% TESTS
function tests


%% KD-TREE FROM A LIST OF POINTS
clearvars,clc

nP = 1000000 ; nCoord = 2 ;
P = rand(nP,nCoord) ;
div = 10 ;
order = max(floor(log(nP)/log(div))-1,0) ; 4 ;

profile on
tic
tree = pkg.graph.tree.KD(P,div,order) ;
toc
profile off

%%
clf ;
axis equal
IDX = {tree.Nodes(tree.endNodes).Data} ;
colors = rand(numel(IDX),3) ; get(gca,'colororder') ;
for ii = 1:numel(IDX)
    cc = mod(ii-1,size(colors,1))+1 ;
    patch('vertices',P(IDX{ii},:),'faces',(1:numel(IDX{ii}))','facecolor','none','marker','.','markersize',5,'markeredgecolor',colors(cc,:)) ;
end


%% TEST POINTS IN BOXES
clearvars
%clc

nCoord = 2 ;
nP = 1000000 ;
nB = 10000 ; szB = 0.1*(1/nB)^(1/nCoord) ; tol = 0.00 ;
div = 10 ;
order = min(ceil(log(nB)/log(div))-1,12) ;

P = rand(nP,nCoord) ;
B = rand(nB,1,nCoord) + (szB/2*[-1 1]).*rand(nB,1,nCoord) ;

boxes = permute(B,[1 3 2]) ; % [nB nCoord 2] ;
bC = mean(boxes,3) ;

encBox = @(ii)cat(3,min(boxes(ii,:,1),[],1)-tol,max(boxes(ii,:,2),[],1)+tol) ;
inBox = @(p,box)all(p>=box(:,:,1),2) & all(p<=box(:,:,2),2) ;

profile on

tic
tree = pkg.graph.tree.KD(bC,div,order) ;
boxTree = tree.nodeFcn(encBox) ; 
inBoxTree = boxTree.applyTest(@(ii,box)inBox(P(ii,:),box),nP) ;
%
toc

profile off
%%
nPtsInNode = cellfun(@numel,{inBoxTree.Nodes.Data})' ;
nPtsInLayer = accumarray(tree.nodeLayer',nPtsInNode') ;
nTests = sum(nPtsInNode(tree.endNodes).*cellfun(@numel,{tree.Nodes(tree.endNodes).Data})')
ratio = nP*nB/nTests

% run final tests
ttt = tic ;
endNodes = tree.endNodes ;
% box indices
    ib = {tree.Nodes(endNodes).Data} ;
    nb = cellfun(@numel,ib) ;
% point indices
    ip = {inBoxTree.Nodes(endNodes).Data} ;
    np = cellfun(@numel,ip) ;
% index repetition
    ip = repelem(cat(2,ip{:}),repelem(nb,np)) ;
    ib = repelem(ib,np) ; ib = cat(2,ib{:}) ;
% test
    in = inBox(P(ip,:),boxes(ib,:,:)) ;
    ip = ip(in) ; ib = ib(in) ;
toc(ttt)

profile off

%%
clf ;
axis equal
patch('vertices',P,'faces',(1:nP)','facecolor','none','marker','.','markeredgecolor','r')
bP = @(box) [box(:,:,1) ; [box(:,1,2) box(:,2,1)] ; box(:,:,2) ; [box(:,1,1) box(:,2,2)]] ;
patch('vertices',bP(boxes),'faces',[(1:nB)' + (0:3)*nB],'facecolor','none','edgecolor','k') ;
layer = boxTree.nodeLayer ;
colors = hsv(max(layer)) ;
for ll = 1:max(layer)
    bbb = cat(1,boxTree.Nodes(layer==ll).Data) ;
    patch('vertices',bP(bbb),'faces',[(1:size(bbb,1))' + (0:3)*size(bbb,1)],'facecolor','none','edgecolor',colors(ll,:),'linewidth',1.5,'linestyle','-.') ;
end



end

