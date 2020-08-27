classdef Node < handle
%NODE A tree node
    
properties
    Data % included data in the node
end
properties (AbortSet)
    Tree = pkg.data.tree.Node.empty
    Parent = pkg.data.tree.Node.empty % parent node in the tree
    Children = pkg.data.tree.Node.empty % children nodes in the tree
end

methods
    function this = Node(tree,parent,data)
    % Constructor
        if nargin>=1 ; this.Tree = tree ; end
        if nargin>=2 ; this.Parent = parent ; end
        if nargin>=3 ; this.Data = data ; end
    end
    
    function delete(this)
    % Destructor
        if isvalid(this.Tree)
            this.Tree.Nodes(this.Tree.Nodes==this) = [] ;
            if ~isempty(this.Parent)
                this.Parent.Children([this.Parent.Children]==this) = [] ;
            end
            [this.Children.Parent] = deal(this.Parent) ;
        end
    end
    
    function set.Parent(this,newParent)
    % Set the node's parent
        if ~isa(newParent,'pkg.data.tree.Node') ; error('Parent must be a pkg.data.tree.Node') ; end
        if isempty(newParent) ; newParent = pkg.data.tree.Node.empty ; return ; end
        if ~isempty(this.Parent)
            this.Parent.Children([this.Parent.Children]==this) = [] ;
        end
        this.Parent = newParent ;
        this.Parent.Children(end+1) = this ;
    end
end

end

