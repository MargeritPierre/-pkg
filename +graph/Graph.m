classdef Graph
%GRAPH Representation of a (directed) graph with nodes and edges
    
%% PROPERTIES
properties
    Edges(:,2) uint32
    Nodes(1,:) pkg.graph.Node
    Directed(1,1) logical = true
end

%% INFOS
methods
    function n = nEdges(this) ; n = size(this.Edges,1) ; end
    function n = nNodes(this) ; n = max(double(max([this.Edges(:) ; 0])),numel(this.Nodes)) ; end
end

%% CONSTRUCTORS
methods
    function this = Graph(varargin)
    % Build the graph
        if mod(nargin,2) ; error('Wrong number of arguments') ; end
        for ii=1:2:nargin
            this.(varargin{ii}) = varargin{ii+1} ;
        end
    end
end

%% SET/GET INTERFACES
methods
    function this = set.Edges(this,edg)
    % Modify the list of edges
        this.Edges = edg ;
        % <TODO> modify nodes ?
    end
end

%% MATRIX REPRESENTATION
methods
    function M = incMat(this)
    % Return the sparse, signed incidence matrix
    % Signed if the graph is directed
    % Mij = -1 if jth edge starts with ith node
    % Mij = 1 if jth edge ends with ith node
        nodIdx = double(this.Edges) ;
        edgIdx = repmat((1:this.nEdges)',[1 2]) ;
        values = repmat([-1 1],[this.nEdges 1]) ;
        if ~this.Directed ; values = abs(values) ; end
        M = sparse(nodIdx(:),edgIdx(:),values(:),this.nNodes,this.nEdges) ;
    end
    
    function M = adjMat(this)
    % Return the sparse, signed adjacency matrix
    % such that Mij = true if one edge is [j i]
    % undirected graph: M is symmetric
        edges = double(this.Edges) ;
        if ~this.Directed ; edges = [edges ; flip(edges,2)] ; end
        M = sparse(edges(:,2),edges(:,1),1,this.nNodes,this.nNodes) ;
    end
    
    function M = degMat(this,dir)
    % Return the sparse, diagonal degree matrix
    % Mii = degree of node (number of attached edges)
    % dir: directions of the edges taken into account
    % dir = 'both' (undirected graph), 'in' (indegree) or 'out' (outdegree)
        if nargin<2 ; dir = defaultDegreeDirection(this) ; end
        deg = nodeDegree(this,dir) ;
        M = spdiags(deg(:),0,this.nNodes,this.nNodes) ;
    end
    
    function dir = defaultDegreeDirection(this)
    % Return the default edge degree counting direction
        if this.Directed ; dir = 'in' ;
        else ; dir = 'both' ;
        end
    end
    
    function M = lapMat(this,dir)
    % Return the sparse Laplacian matrix L = D(dir)-A
    % where D(dir) is the degree matrix
    % and A is the adjacency matrix
        if nargin<2 ; dir = defaultDegreeDirection(this) ; end
        M = this.degMat(dir)-this.adjMat ;
    end
end

%% NODE FEATURES
methods
    function deg = nodeDegree(this,dir)
    % Return the node degree corresponding to a given direction
        if nargin<2 ; dir = defaultDegreeDirection(this) ; end
        N = this.nNodes ;
        switch dir
            case 'both'
                idx = this.Edges(:) ;
            case 'in'
                idx = this.Edges(:,2) ;
            case 'out'
                idx = this.Edges(:,1) ;
        end
        deg = accumarray(idx,1,[N 1],[],[],true) ;
    end
    
    function nn = startNodes(this)
    % Return start nodes (nodes with no parents)
        nn = ~ismember(1:this.nNodes,this.Edges(:,2)') ;
    end

    function nn = endNodes(this)
    % Return end nodes (nodes with no children)
        nn = ~ismember(1:this.nNodes,this.Edges(:,1)') ;
    end

    function nn = loneNodes(this)
    % Return isolated nodes (nodes with no associated edges)
        nn = ~ismember(1:this.nNodes,this.Edges(:)') ;
    end
end
   
%% WALK IN THE GRAPH
methods
    function nodIdx = neightborNodes(this,dir)
    % Return the indices of nehgtboring nodes corresponding to a given direction
    % nodIdx is a cell array of size [nNodes 1]
        if nargin<2 ; dir = defaultDegreeDirection(this) ; end
    % Edge indices to use
        switch dir
            case 'both' % all attached nodes
                starts = this.Edges(:) ;
                ends = [this.Edges(:,2) ; this.Edges(:,1)] ;
            case 'in' % previous nodes
                starts = this.Edges(:,2) ;
                ends = this.Edges(:,1) ;
            case 'out' % next nodes
                starts = this.Edges(:,1) ;
                ends = this.Edges(:,2) ;
        end
    % Sort the indices
        [starts,is] = sort(starts,'ascend') ;
        ends = ends(is) ;
    % Put in a cell
        nNodIdx = accumarray(starts(:),1,[this.nNodes 1]) ;
        nodIdx = mat2cell(ends,nNodIdx) ;
    end
    
    function prevIdx = previousNodes(this,nn)
    % Return the nodes attached to the given nodes by an INgoing edge
    % asCell (bool)
        if this.Directed ; dir = 'in' ; else ; dir = 'both' ; end
        prevIdx = this.neightborNodes(dir) ;
        if nargin>1 ; prevIdx = prevIdx(nn) ; end
    end
    
    function nextIdx = nextNodes(this,nn)
    % Return the nodes attached to the given nodes by an OUTgoing edge
        if this.Directed ; dir = 'out' ; else ; dir = 'both' ; end
        nextIdx = this.neightborNodes(dir) ;
        if nargin>1 ; nextIdx = nextIdx(nn) ; end
    end
    
    function dist = distanceTo(this,nn)
    % Return the shortest distance from each node to a list of starting nodes
    % nn can be logical or a list of indices
        if nargin<2 ; nn = this.startNodes ; end
        if islogical(nn) ; nn = find(nn) ; end
        dist = NaN(1,this.nNodes) ;
        dist(this.loneNodes) = Inf ;
        dist(nn) = 0 ;
        dd = 0 ;
        while any(isnan(dist))
            dd = dd+1 ;
            previousNodes = find(dist==dd-1) ;
            availableNodes = find(isnan(dist)) ;
            nextEdg = ismember(this.Edges(:,1)',previousNodes) ;
            nextNodes = ismember(availableNodes,this.Edges(nextEdg,2)') ;
            dist(availableNodes(nextNodes)) = dd ;
        end
    end
end

%% GRAPH TOPOLOGY
methods
    function [b,nodeIdx,order] = isCyclic(this)
    % Return true if the graph is cyclic
    % uses consecutive powers of the adjacency matrix (!!!)
        if ~this.Directed % an undirected graph is obviously cyclic
            b = true ; nodeIdx = 1:this.nNodes ; order = 2 ;
            return
        end
        A = this.adjMat ;
        An = A ;
        nodeIdx = [] ; order = 1 ;
        maxOrder = min(this.nEdges,this.nNodes) ;
        while order<maxOrder && ~any(nodeIdx)
            order = order + 1 ;
            An = An*A ;
            nodeIdx = logical(diag(An)) ;
        end
        nodeIdx = find(nodeIdx) ;
        b = ~isempty(nodeIdx) ;
        if ~b ; order = Inf ; end
    end
    
    function C = connComp(this)
    % Return the connected components of the graph
    % -> nullspace of the laplacian matrix
    % C: [nNodes nComp] logical
        dir = defaultDegreeDirection(this) ;
        L = full(this.lapMat(dir)) ;
        C = null(L) ;
    end
    
    function sp = sparsestCut(this)
    % Return the sparsest cut logical indices
    % sp is logical [nNodes 1] (gives the side of the cut)
        dir = defaultDegreeDirection(this) ;
        L = this.lapMat(dir) ;
        [W,~] = eigs(L,2,'sm') ;
        sp = W(:,2) ;
    end
end

%% EULERIAN PATHS
% paths that go trought every edge once and ONLY once 
% (vertices can be attained more than once)
methods
    function bool = isEulerian(this)
    % Is the graph Eulerian ?
        if this.Directed % directed graphs
        % at most 2 nodes can have different in & out degrees
            inDeg = this.nodeDegree('in') ;
            outDeg = this.nodeDegree('out') ;
            diffDeg = inDeg-outDeg ;
            bool = sum(abs(diffDeg))<=2 && sum(sign(diffDeg))==0 ;
        else % undirected graphs
        % at most 2 nodes can have uneven degreee
            deg = this.nodeDegree ;
            iseven = mod(deg,2)==0 ;
            bool = sum(iseven)<=2 ;
        end
    end
end

%% OPERATIONS ON GRAPH
methods
    function G = subGraph(this,nn)
    % Extract a subpart of the graph with given nodes nn
    % nn can be either logical or indicial
        if islogical(nn) ; nn = find(nn) ; end
        G = this ; % Create a copy
        if ~isempty(G.Nodes) % Copy node data if needed
            G.Nodes = G.Nodes(min(nn,numel(G.Nodes))) ;
        end
    % Cull unused edges and change node indices
        [~,newIdx] = ismember(this.Edges,nn) ;
        keepEdg = all(newIdx,2) ;
        G.Edges = newIdx(keepEdg,:) ;
    end
end

%% GRAPHICAL REPRESENTATION
methods
    function X = autoNodeCoordinates(this)
    % Return a set of node coordinates for plotting
        if 0 % random distribution
            X = rand(this.nNodes,2) ;
        else % circle
            t = (0:this.nNodes-1)'/this.nNodes*2*pi ;
            X = [cos(t) sin(t)] ;
        end
    end
    
    function h = plot(this,varargin)
    % Plot the tree in the current axes
        h = pkg.graph.GraphPlot(varargin{:}) ;
        h.Parent = gca ;
        h.Graph = this ;
    end
end

end





%% TEST FUNCTIONS
function tests

%% SPEED TESTS
clearvars
clc
nNodes = 10000 ; nEdges = 10000 ;
edges = unique(randi(nNodes,[nEdges 2]),'rows') ;

profile on
tic ;
G = pkg.graph.Graph('Directed',false,'Edges',edges) ;
sn = G.startNodes ;
en = G.endNodes ;
ln = G.loneNodes ;
idx = G.neightborNodes ;
I = G.incMat ;
A = G.adjMat ;
D = G.degMat ;
L = G.lapMat ;
toc ;
profile off

% clf
% pl = plot(G) ;
% pl.ShowLabels = {'Nodes' 'Edges'} ;
% axis equal


%% CONNECTED COMPONENTS
clearvars
clc

topo = 'rand' ;
switch topo
    case 'rand'
        nNodes = 10 ; nEdges = 10 ;
        edges = unique(randi(nNodes,[nEdges 2]),'rows') ;
        edges(range(edges,2)==0,:) = [] ;
    case 'loop'
        N = 9 ;
        edges = [1:N ; circshift(1:N,1)]' ;
    case 'loops'
        nLoops = 5 ; nodeRange = [3 10] ;
        edges = [] ; nNodes = 0 ;
        for ll = 1:nLoops
            N = randi(nodeRange) ;
            edges = [edges ; nNodes + [1:N ; circshift(1:N,-1)]' ] ;
            nNodes = nNodes + N ;
        end
    case '8'
        nLoops = 3 ; nodeRange = [10 10] ;
        edges = [] ; nNodes = 0 ;
        N = randi(nodeRange,[1 nLoops]) ;
        for ll = 1:nLoops
            edges = [edges ; nNodes + [1:N(ll) ; circshift(1:N(ll),-1)]' ] ;
            nNodes = nNodes + N(ll) ;
        end
        firstNodes = [0 cumsum(N(1:end-1))]+1 ;
        if 0 % weld nodes
            edges(firstNodes,:) = 1 ;
        else % add an edge
            edges = [edges ; [firstNodes(1:end-1) ; firstNodes(2:end)]'] ;
        end
    case 'custom'
        edges = [1 2 ; 2 3 ; 3 1 ; 1 4 ; 4 5 ; 5 6] ; % spoon
        %edges = [1 2 ; 2 3 ; 3 4] ; % line
end

G = pkg.graph.Graph('Directed',true,'Edges',edges) ;
tic
[b,nodeIdx,order] = G.isCyclic
C = G.connComp
e = G.isEulerian
%sp = G.sparsestCut
toc

clf
pl = plot(G) ;
pl.ShowLabels = {'Nodes' 'Edges'} ;
axis equal

%%
end

