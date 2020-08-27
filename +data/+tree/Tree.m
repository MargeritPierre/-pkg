classdef Tree < handle
%TREE A data tree
    
properties
    Nodes = pkg.data.tree.Node.empty
end

methods
    function this = Tree(varargin)
    % Constructor
    end
    
    function delete(this)
    % Destructor
        delete(this.Nodes) ;
    end
    
    function rn = rootNodes(this)
    % Return root node indices (nodes with no parent node)
        rn = cellfun(@isempty,{this.Nodes.Parent}) ;
    end
    
    function en = endNodes(this)
    % Return end node indices (nodes with no children nodes)
        en = cellfun(@isempty,{this.Nodes.Children}) ;
    end
    
    function n = nNodes(this)
    % Number of nodes in the tree
        n = numel(this.Nodes) ;
    end
    
    function lay = deepness(this)
    % Deepness of each node
        ll = 0 ;
        lay = zeros(1,this.nNodes) ;
        indLay = this.rootNodes ;
    end
end

end

function tests
%% CREATE A TREE WITH RANDOM NODES connectivities
clearvars
N = 1000 ;
tree = pkg.data.tree.Tree ;

profile on
tic ;
for nn = 1:N
    newNode = pkg.data.tree.Node(tree) ;
    if numel(tree.Nodes)>0
        newNode.Parent = tree.Nodes(randi(numel(tree.Nodes))) ;
    end
    tree.Nodes(end+1) = newNode ;
end
toc
profile off

%%
end

