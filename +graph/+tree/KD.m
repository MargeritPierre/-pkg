function tree = KD(data,order)
%KD Build a kd tree from data (median-splitted branches)
% data: [nData nCoord nFeatures]
% examples: 
%   - nd-values: [nValues nCoord 1]
%   - nd-domains: [nDomains nCoord 2] cat(3,lower_nd-values,upper_nd-values)

if nargin==0 ; error('Not enough arguments') ; end

% Infos
data = data(:,:,:) ; % [nData nCoord nFeatures]
[nData,nCoord,nFeatures] = size(data) ;
permData = permute(data,[1 3 2]) ; % [nData nFeatures nCoord]

% Automatic order: maximum order tree
if nargin<2 ; order = floor(log2(nData)) ; end

% Init tree
tree = pkg.data.tree.Tree ;

% Root node with all data indices
tree.Nodes = pkg.data.tree.Node(tree) ;
tree.Nodes.Data = 1:nData ;

% Global bbox
allValues = reshape(permData,[nData*nFeatures nCoord]) ;
bbox = [min(allValues,[],1) ; max(allValues,[],1)] ;

% Loop
prevLayerNodes = tree.Nodes ;
for oo = 1:order
    currentNodes = pkg.data.tree.Node.empty ;
    coord = mod(oo-1,nCoord)+1 ;
    for nn = 1:numel(prevLayerNodes)
        node = prevLayerNodes(nn) ;
        ind = node.Data ;
        nodeData = reshape(data(ind,coord,:),[],nFeatures) ; % [nDataInNode nFeatures]
        med = median(nodeData,2) ;
        higher = med>=median(med) ;
        currentNodes(end+1) = pkg.data.tree.Node(tree,node,ind(~higher)) ;
        currentNodes(end+1) = pkg.data.tree.Node(tree,node,ind(higher)) ;
    end
    tree.Nodes = [tree.Nodes currentNodes] ;
    prevLayerNodes = currentNodes ;
end

end

%% TESTS
function tests
%% KD-TREE FROM A LIST OF POINTS
clearvars
P = rand(10000,2) ;
profile on
tic
tree = pkg.data.tree.KD(P) ;
toc
profile off


end

