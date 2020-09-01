classdef Tree < pkg.graph.Graph
%TREE A data tree
    
properties
end

methods
    function this = Tree(varargin)
    % Constructor
        this@pkg.graph.Graph(varargin{:}) ;
    end
    
    function lvl = nodeLayer(this)
    % Return indices corresponding to the node level (or layer)
        lvl = this.distanceTo(this.startNodes)+1 ;
    end
end

end

