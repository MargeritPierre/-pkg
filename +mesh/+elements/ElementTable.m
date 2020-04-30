classdef ElementTable
%ELEMENTTABLE Element information for meshes
% Works for elements, faces and edges
    
%% ELEMENT INFORMATIONS
properties (AbortSet)
<<<<<<< HEAD
    % Type of elements [nElems 1]
=======
    % Unique list of element types
>>>>>>> parent of 420de96... New Commit
    Types pkg.mesh.elements.AbstractElement
    % Element-node connectivities [nElems nMaxNodeByElem]
    % Zeros denote invalid indices
<<<<<<< HEAD
    NodeIndices uint32
=======
    % Format : [elmtTypIdx NodesIndices]
    Indices uint32
>>>>>>> parent of 420de96... New Commit
end

%% TABLE INFORMATIONS
methods  
    % Number of elements in the table
    function val = nElems(this) ;  val = size(this.NodeIndices,1) ; end
    % Number of nodes in each element
<<<<<<< HEAD
    function val = nNodes(this) ; val = sum(this.NodeIndices>0,2) ; end
    % Maximum number of nodes by element
    function val = nMaxNodesByElem(this) ;  val = size(this.NodeIndices,2) ; end
    % Unique list of node indices
    function val = uniqueNodeIndices(this) ; val = unique(this.NodeIndices(this.NodeIndices~=0)) ; end
=======
    function val = nNodes(this) ; val = sum(this.Indices>0,2) ; end
    % Maximum number of nodes by element
    function val = nMaxNodesByElem(this) ;  [~,val] = cellfun(@size,{this.Indices}) ; end
    % Unique list of node indices
    function val = uniqueIndices(this) ; val = unique(this.Indices(this.Indices~=0)) ; end
>>>>>>> parent of 420de96... New Commit
    % Unique list of element types
    function [types,elmtIdx] = uniqueTypes(this) ; [types,~,elmtIdx] = unique(this.Types) ; end
end

%% CONSTRUCTOR
methods
    function this = ElementTable(varargin)
    % Class constructor
    % Process arguments
        % Odd number of arguments ?
        if mod(nargin,2)==1 ; error('Wrong number of arguments') ; end
        if nargin>0
            PropNames = varargin(1:2:end-1) ; 
            PropValues = varargin(2:2:end) ; 
            for pp = 1:numel(PropNames)
                this.(PropNames{pp}) = PropValues{pp} ;
            end
        end
    end
end

%% SET / GET interface
methods
    function this = set.Types(this,types)
    % Set the Types property
        this.Types = types(:) ;
    % Check data consistency
        if numel(types)>this.nElems
        % Too many types, add empty elements
<<<<<<< HEAD
            this.NodeIndices(end+1:numel(types),:) = 0 ;
        elseif numel(types)<this.nElems 
        % Not enough types, try to find each type by node number
            this.Types = this.findElementTypes ;
        end
    end
    
    function this = set.NodeIndices(this,indices)
    % Set the NodeIndices property
        this.NodeIndices = indices(:,:) ;
    % Check data consistency
        if numel(this.Types)==0
        % No element type has been provided, fill with abstract elements
=======
            this.Indices(end+1:numel(types),:) = 0 ;
        elseif numel(types)<this.nElems 
        % Not enough types, try to find each type by node number
            this = this.assignElementTypes ;
        end
    end
    
    function this = set.Indices(this,indices)
    % Set the Indices property
        this.Indices = indices(:,:) ;
    % Check data consistency
        if numel(this.Types)==0
        % No element type has been provided, fill with abstract elements
            %this.Types(1:this.nElems) = pkg.mesh.elements.EmptyElement ;
>>>>>>> parent of 420de96... New Commit
            this.Types = repmat(pkg.mesh.elements.EmptyElement,[this.nElems 1]) ;
        elseif this.nElems>numel(this.Types)
        % Not enough element types are provided
        % assign auto element types on missing lines only !
            types = this.findElementTypes ;
            indMissing = numel(this.Types)+1:this.nElems ;
            this.Types(indMissing) = types(indMissing) ;
        elseif this.nElems<numel(this.Types) 
        % Too many elements are present, remove the ones that are not used
            this.Types = this.Types(1:this.nElems) ;
        end
<<<<<<< HEAD
    end
    
    function types = findElementTypes(this)
    % Return the list of element types corresponding to the current list of
    % NodeIndices. Use the number of nodes by element. 
    % Unique element type list
        types = this.uniqueTypes ; 
    % Number of nodes in each element type
        nNodesInElements = [types.nNodes] ;
        if numel(nNodesInElements)~=numel(unique(nNodesInElements))
            warning('Some elements types have the same node number: automatic type assignment is ambiguous.') ;
        end
    % Number of nodes in each line of the node indices
        nNodes = this.nNodes ;
    % Initialize with empty elements by default
        newTypes = repmat(pkg.mesh.elements.EmptyElement,[this.nElems 1]) ;
     % Assign with the count the number of nodes
        for tt = 1:numel(types)
            newTypes(nNodesInElements(tt)==nNodes) = types(tt) ;
        end
        types = newTypes ;
=======
    end
    
    function types = findElementTypes(this)
    % Return the list of element types corresponding to the current list of
    % Indices. Use the number of nodes by element. 
    % Unique element type list
        types = this.uniqueTypes ; 
    % Number of nodes in each element type
        nNodesInElements = [types.nNodes] ;
        if numel(nNodesInElements)~=numel(unique(nNodesInElements))
            warning('Some elements types have the same node number: automatic type assignment is ambiguous.') ;
        end
    % Number of nodes in each line of the node indices
        nNodes = this.nNodes ;
    % Initialize with empty elements by default
        newTypes = repmat(pkg.mesh.elements.EmptyElement,[this.nElems 1]) ;
     % Assign with the count the number of nodes
        for tt = 1:numel(types)
            newTypes(nNodesInElements(tt)==nNodes) = types(tt) ;
        end
        types = newTypes ;
    end
end


%% CONNECTIVITY RETRIEVING
methods
    function [edges,ie] = uniqueEdges(this)
    % Return the list of unique edges
    end
    
    function edges = allEdges(this)
    % Return the complete list of edges (with duplicates)
    % Build Edge Table
        edges = [this.Types.Edges] ; % Local indices (ELEMENT node number)
    % Globalize indices
        nEdges = [this.Types.nEdges] ; % Edges by element
        edgElem = repelem(1:this.nElems,nEdges) ; % Element associated to each edge
        edgElem = repmat(edgElem(:),edges.nMaxNodesByElem) ;
        valid = edges.Indices>0 ; % Valid edge indices
        iii = sub2ind(size(this.Indices),edgElem(valid),edges.Indices(valid)) ;
        edges.Indices(valid) = this.Indices(iii) ; % Global indices (MESH node number)
    end
end


%% INDICES MANIPULATION
methods
    function indices = paddedIndices(this,N,indices)
    % Pad the table indices with zeros to match a given N
        if nargin<3 ; indices = this.Indices ; end
        N0 = size(indices,2) ;
        if N0~=N ; indices = [indices(:,1:min(N0,N)) zeros(this.nElems,N-N0,'uint32')] ; end
    end
    
    function indices = catIndices(tables)
    % Concatenate indices of multiple tables
    % Retrieve indices
        IND = {tables.Indices} ;
    % Retrieve indices size
        [nElems,nMax] = cellfun(@size,IND) ;
        [uNMax,~,ia] = unique(nMax) ;
    % If all lengths match
        if numel(uNMax)==1
            indices = cat(1,IND{:}) ; 
            return ;
        end
    % Initialize
        indices = zeros(sum(nElems),max(nMax),'uint32') ; 
    % Fill the table
        emax = cumsum(nElems) ;
        emin = emax-nElems+1 ;
        ee = arrayfun(@colon,emin,emax,'UniformOutput',false) ;
        for nn = 1:numel(uNMax)
            iii = ia==nn ;
            eee = cat(2,ee{iii}) ;
            indices(eee,1:uNMax(nn)) = cat(1,IND{iii}) ;
        end
    end
    
    function [table,ia] = unique(this)
    % Return a table of UNIQUE elements
        [~,ie,ia] = unique(sort(this.Indices,2),'rows','stable') ;
        table = pkg.mesh.elements.ElementTable('Types',this.Types(ie),'Indices',this.Indices(ie,:)) ;
>>>>>>> parent of 420de96... New Commit
    end
end

%% INDEXATION, CONCATENATION, ETC
methods
    
<<<<<<< HEAD
    function varargout = subsref(this,s)
    % Sub-referencing
        switch s(1).type
            case '()'
                s(1).subs = {s(1).subs{1},':'} ; 
                varargout{1} = pkg.mesh.elements.ElementTable(...
                                    'Types',subsref(this.Types,s(1)) ...
                                    ,'NodeIndices',subsref(this.NodeIndices,s(1)) ...
                                    ) ;
                if numel(s)>1
                    [varargout{1:nargout}] = subsref(varargout{1},s(2:end)) ;
                end
            otherwise
                [varargout{1:nargout}] = builtin('subsref',this,s);
        end
=======
    function table = cat(~,varargin)
    % Concatenate element tables to an unique table
    % Convert to tables
        tables = builtin('cat',1,varargin{:}) ; %cellfun(@pkg.mesh.elements.ElementTable,varargin) ;
    % Concatenate Element Types
        types = cat(1,tables.Types) ;
    % Concatenate Indices
        indices = catIndices(tables) ;
    % Create the table
        %table = pkg.mesh.elements.ElementTable('Indices',indices,'Types',types) ; 
        table = pkg.mesh.elements.ElementTable() ; 
        table.Types = types ;
        table.Indices = indices ; 
>>>>>>> parent of 420de96... New Commit
    end
    
    function newTable = vertcat(this,otherTable)
    % Concatenate two tables
        N = max(this.nMaxNodesByElem,otherTable.nMaxNodesByElem) ;
        newNI = [this.NodeIndices zeros(this.nElems,N-this.nMaxNodesByElem,'uint32') ; ...
                 otherTable.NodeIndices zeros(otherTable.nElems,N-otherTable.nMaxNodesByElem,'uint32')] ;
        newTable = pkg.mesh.elements.ElementTable(...
                        'Types',[this.Types(:) ; otherTable.Types(:)] ...
                        ,'NodeIndices',newNI ...
                        ) ;
    end
    
    function newTable = horzcat(this,otherTable) 
    % Concatenate two tables
        newTable = vertcat(this,otherTable) ;
    end
    
end

end
