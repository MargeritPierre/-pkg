classdef ElementTable
%ELEMENTTABLE Element information for meshes
% Works for elements, faces and edges
% Contains all connectivities of the mesh
    
%% ELEMENT INFORMATIONS
properties
    % Unique list of element types
    Types pkg.mesh.elements.AbstractElement
    % Element-node connectivities [nElems nMaxNodeByElem+1]
    % Zeros denote invalid indices
    % Format : [elmtTypeIdx NodesIndices]
    Indices uint32
end
properties (Dependent)
    TypeIdx uint32 % Element type indices (wrt this.Types)
    NodeIdx uint32 % Element Node indices (wrt mesh nodes or element nodes)
end
methods
end

%% TABLE INFORMATIONS
methods  
    % Number of elements in the table
    function val = nElems(this) ;  [val,~] = cellfun(@size,{this.Indices}) ; end
    % Number of nodes in each element
    function val = nNodes(this) ; val = sum(this.NodeIdx>0,2) ; end
    % Maximum number of nodes by element
    function val = nMaxNodesByElem(this) ;  [~,val] = cellfun(@size,{this.Indices}) ; val = val-1 ; end % Minus one because of the presence of TypeIdx
    % Unique list of node indices
    function val = uniqueNodeIdx(this) ; val = unique(this.NodeIdx(this.NodeIdx~=0)) ; end
    % Unique list of element types
    %function [types,elmtIdx] = uniqueTypes(this) ; [types,~,elmtIdx] = unique(this.Types) ; end
end

%% CONSTRUCTOR
methods
    function this = ElementTable(varargin)
    % Class constructor (handle conversion to the type)
    % Process arguments
        if nargin==1 % only one argument (indices)
            data = varargin{1} ;
            switch class(data)
                case 'pkg.mesh.elements.ElementTable'
                    this = data ;
                otherwise
                    this = pkg.mesh.elements.ElementTable('Indices',varargin{1}) ;
            end
        elseif mod(nargin,2)==1 % Odd number of arguments
            error('Wrong number of arguments') ;
        elseif nargin>0
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
    function idx = get.TypeIdx(this)
    % Get element type indices
        idx = this.NodeIdx(:,1) ;
    end
    
    function idx = get.NodeIdx(this)
    % Get element node indices
        idx = this.NodeIdx(:,2:end) ;
    end
    
    function this = set.TypeIdx(this,idx)
    % Set element type indices
        if numel(idx)~=this.nElems ; error('The provided indices must be of size [nElems 1]') ; end
        if any(idx>numel(this.Types)) ; error('The indices must range in the number of available element types') ; end 
        this.Indices = [idx(:) this.NodeIdx] ;
    end
    
    function this = set.NodeIdx(this,idx)
    % Set element node indices
        if size(idx,1)~=this.nElems ; error('The provided list of indices must have [nElems] rows') ; end
        this.Indices = [this.TypeIdx idx(:,:)] ;
    end
    
    function this = set.Types(this,types)
    % Set the list of unique elements
    % Unique list
        [this.Types,~,typeIdx] = unique(types(:)) ;
    % Change the type index list
        if numel(typeIdx)==this.nElems % One type by element, already ok
            this.TypeIdx = typeIdx ;
        else % We have to find which type fit with the given element node indices
            this.TypeIdx = this.assignElementTypeIdx() ;
        end
    end
    
    function this = set.Indices(this,indices)
    % Set the Indices property
    % Check if element type indices are valid
        % In the right range
            valid = indices(:,1)>0 & indices(:,1)<numel(this.Types) ;
        % The number of node indices is correct
            nNodes = reshape([this.Types.nNodes],[],1) ;
            valid(valid) = valid(valid) & sum(indices(valid,2:end)>0,2)==nNodes(indices(valid,1)) ;
    % Auto assign type indices if needed
        if any(~valid) 
            indices(~valid,1) = this.assignElementTypeIdx(indices(~valid,2:end)) ;
        end
    % Set
        this.Indices = indices(:,:) ;
    end
end


%% CONNECTIVITY RETRIEVING
methods
%     function [edges,ie] = uniqueEdges(this)
%     % Return the list of unique edges
%     end
    
%     function edges = allEdges(this)
%     % Return the complete list of edges (with duplicates)
%     % Build Edge Table
%         edges = [this.Types.Edges] ; % Local indices (ELEMENT node number)
%     % Globalize indices
%         nEdges = [this.Types.nEdges] ; % Edges by element
%         edgElem = repelem(1:this.nElems,nEdges) ; % Element associated to each edge
%         edgElem = repmat(edgElem(:),edges.nMaxNodesByElem) ;
%         valid = edges.Indices>0 ; % Valid edge indices
%         iii = sub2ind(size(this.Indices),edgElem(valid),edges.Indices(valid)) ;
%         edges.Indices(valid) = this.Indices(iii) ; % Global indices (MESH node number)
%     end
end


%% INDICES MANIPULATION
methods
    function typeIdx = assignElementTypeIdx(this,nodeIdx,types)
    % Return the list of element types indices corresponding to a list of Node indices. 
    % Uses the number of nodes by element. 
        if nargin<2 ; nodeIdx = this.NodeIdx ; end
        if nargin<3 ; types = this.Types ; end
    % Number of nodes in each available element type
        nNodeInType = [types.nNodes] ;
        if numel(nNodeInType)~=numel(unique(nNodeInType))
            warning('Some elements types have the same node number: automatic type assignment is ambiguous.') ;
        end
    % Find new TypeIdx
        nValidNodeIdx = sum(nodeIdx>0,2) ;
        [~,typeIdx] = ismember(nValidNodeIdx,nNodeInType) ;
    end
    
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
        table = pkg.mesh.elements.ElementTable('Types',this.Types,'Indices',this.Indices(ie,:)) ;
    end
end


%% INDEXATION, CONCATENATION, ETC
methods
%     function varargout = subsref(this,s)
%     % Sub-referencing
%     % Allow to extract a sub-part of ONE table with element indices
%         if numel(this)>1 
%             [varargout{1:nargout}] = builtin('subsref',this,s) ; 
%             return ;
%         end
%     % Sub-table extraction
%         switch s(1).type
%             case '()'
%                 if this.nElems==0 % No elements to extract
%                     varargout{1} = pkg.mesh.elements.ElementTable() ;
%                 else % There is elements to extract
%                     s(1).subs = {s(1).subs{1},':'} ; 
%                     varargout{1} = pkg.mesh.elements.ElementTable(...
%                                         'Types',subsref(this.Types,s(1)) ...
%                                         ,'Indices',subsref(this.Indices,s(1)) ...
%                                         ) ;
%                 end
%                 if numel(s)>1
%                     [varargout{1:nargout}] = subsref(varargout{1},s(2:end)) ;
%                 end
%             otherwise
%                 [varargout{1:nargout}] = builtin('subsref',this,s);
%         end
%     end
    
    function table = cat(~,varargin)
    % Concatenate element tables to an unique table
    % Convert to tables
        tables = builtin('cat',1,varargin{:}) ; %cellfun(@pkg.mesh.elements.ElementTable,varargin) ;
    % Concatenate Element Types
        types = {tables.Types} ;
        nTypes = cellfun(@numel,types) ;
    % Concatenate Indices
        indices = catIndices(tables) ;
    % Create the table
        %table = pkg.mesh.elements.ElementTable('Indices',indices,'Types',types) ; 
        table = pkg.mesh.elements.ElementTable() ; 
        table.Types = types ;
        table.Indices = indices ; 
    end
    
    function table = vertcat(varargin)
        table = cat(1,varargin{:}) ;
    end

    function table = horzcat(varargin)
        table = cat(1,varargin{:}) ;
    end
end


%% CONNECTIVITY LISTS

end
