classdef Tree < pkg.graph.Graph
%TREE A data tree
    
properties
end

%% CONSTRUCTORS
methods
    function this = Tree(varargin)
    % Constructor
        this@pkg.graph.Graph(varargin{:}) ;
    end
end

%% TREE CONNECTIVITY
methods
    function lvl = nodeLayer(this)
    % Return indices corresponding to the node level (or layer)
        lvl = this.distanceTo(this.startNodes)+1 ;
    end
    
    function n = nLayers(this)
    % Return the number of layers in the tree
        n = max(this.nodeLayer) ;
    end
    
    function n = nNodesByLayers(this)
    % Retrun the number of nodes in each layer of the tree
        n = accumarray(this.nodeLayer',1) ;
    end
    
    function parents = parentNode(this,nn)
    % Return the parent node of given node indices
        if nargin<2 ; nn = 1:this.nNodes ; end
        if islogical(nn) ; nn = find(nn) ; end
        [valid,inEdg] = ismember(nn,this.Edges(:,2)) ;
        parents = zeros(1,numel(nn)) ;
        parents(valid) = this.Edges(inEdg(valid),1) ;
    end
    
    function children = childNodes(this,nn)
    % Return the parent node of given node indices
        if nargin<2 ; nn = 1:this.nNodes ; end
        children = this.nextNodes(nn) ;
    end
    
    function rn = rootNode(this)
    % Return the root node of the tree
        rn = find(this.startNodes) ;
    end
end

%% FUNCTION DISPATCHING
methods
    function tree = nodeFcn(this,fcn)
    % Apply a function to the node data
        tree = this ; % copy the tree
        out = cellfun(fcn,{tree.Nodes.Data},'UniformOutput',false) ;
        [tree.Nodes.Data] = deal(out{:}) ;
    end
    
    function tree = applyTest(this,testFcn,nData)
    % Apply a test using the tree data structure
    % testFcn is an anonymous function handle taking two arguments:
    %   - indices returned by a previous test
    %   - the node data
        tree = this ; % make a copy
    % Tree structure
        nodeLayer = this.nodeLayer ;
        parentNodes = this.parentNode ;
    % Indices
        IDX = cell(1,tree.nNodes) ;
        for ll = 1:max(nodeLayer)
            nn = find(nodeLayer==ll) ;
            if ll==1
                ii = repmat({1:nData},[1 numel(nn)]) ;
            else
                ii = IDX(parentNodes(nn)) ;
            end
            data = {this.Nodes(nn).Data} ;
            bool = cellfun(testFcn,ii,data,'UniformOutput',false) ;
            ii = cellfun(@(ii,bool)ii(bool),ii,bool,'UniformOutput',false) ;
            [IDX{nn}] = deal(ii{:}) ;
        end
    % Put indices in the tree
        [tree.Nodes.Data] = deal(IDX{:}) ;
    end
    
    function tree = applyTest_bkp(this,testFcn,nData)
    % Apply a test using the tree data structure
    % testFcn is an anonymous function handle taking two arguments:
    %   - indices returned by a previous test
    %   - the node data
        tree = this ; % make a copy
    % Tree structure
        rootNode = this.rootNode ;
        [~,nodeOrder] = sort(this.nodeLayer,'ascend') ;
        parentNodes = this.parentNode ;
    % Indices
        IDX = cell(tree.nNodes,1) ;
        for nn = nodeOrder
            if nn==rootNode ; ii = 1:nData ; % test vs. all data
            else ; ii = IDX{parentNodes(nn)} ; % test vs. parent node data
            end
            if isempty(ii) ; continue ; end
            data = this.Nodes(nn).Data ;
            bool = testFcn(ii,data) ;
            IDX{nn} = ii(bool) ;
        end
    % Put indices in the tree
        [tree.Nodes.Data] = deal(IDX{:}) ;
    end
end

%% GRAPHICAL REPRESENTATION
methods
    function h = plot(this,ax)
    % Plot the tree in the current axes
        if nargin<2 ; ax = gca ; end
        h = hggroup(ax) ;
        layer = this.nodeLayer ;
        nNodesInLayer = accumarray(layer(:),1) ;
        nNodesInPreviousLayers = [0 cumsum(nNodesInLayer(1:end-1))'] ;
        x = (1:this.nNodes) - nNodesInPreviousLayers(layer) ;% - nNodesInLayer(layer)'/2 ;
        x = x./nNodesInLayer(layer)' - 0.5./nNodesInLayer(layer)' ;
        x = x*max(nNodesInLayer) ;
        X = [x(:) layer(:)] ;
        p = patch('vertices',X,'Faces',this.Edges,'Parent',h) ;
        p.FaceColor = 'none' ;
        p.EdgeColor = 'k' ;
        p.Marker = 'o' ;
        p.MarkerEdgeColor = 'k' ;
        p.MarkerFaceColor = 'w' ;
        p.MarkerSize = 10 ;
    end
end

end

