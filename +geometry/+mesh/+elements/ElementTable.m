classdef ElementTable
%ELEMENTTABLE Element information for meshes
% Works for elements, faces and edges
% Contains all connectivities of the mesh
    
%% ELEMENT INFORMATIONS
properties
    % List of element types
    Types pkg.geometry.mesh.elements.AbstractElement
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
    % Number of elements types in the list
    function val = nTypes(this) ;  val = cellfun(@numel,{this.Types}) ; end
    % Number of nodes in each element
    function val = nNodes(this) ; val = sum(this.NodeIdx>0,2) ; end
    % Number of edges in each element
    function val = nEdges(this) 
        nE = [0 [this.Types.nEdges]] ; 
        val = nE(this.TypeIdx+1)' ; 
    end
    % Number of faces in each element
    function val = nFaces(this) 
        nF = [0 [this.Types.nFaces]] ; 
        val = nF(this.TypeIdx+1)' ; 
    end
    % Maximum number of nodes by element
    function val = nMaxNodesByElem(this) ;  [~,val] = cellfun(@size,{this.Indices}) ; val = val-1 ; end % Minus one because of the presence of TypeIdx
    % Unique list of node indices
    function val = uniqueNodeIdx(this) ; val = unique(this.NodeIdx(this.NodeIdx~=0)) ; end
    % Is the table empty
    function empty = isempty(this) ; empty = builtin('isempty',this) || this.nElems==0 ; end
end

%% CONSTRUCTOR
methods
    function this = ElementTable(varargin)
    % Class constructor (handle conversion to the type)
    % Process arguments
        if nargin==1 % only one argument (indices)
            data = varargin{1} ;
            switch class(data)
                case 'pkg.geometry.mesh.elements.ElementTable'
                    this = data ;
                otherwise
                    this = pkg.geometry.mesh.elements.ElementTable('NodeIdx',varargin{1}) ;
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
        if isempty(this.Indices)
            idx = [] ;
        else
            idx = this.Indices(:,1) ;
        end
    end
    
    function idx = get.NodeIdx(this)
    % Get element node indices
        if isempty(this.Indices)
            idx = [] ;
        else
            idx = this.Indices(:,2:end) ;
        end
    end
    
    function this = set.TypeIdx(this,idx)
    % Set element type indices cannot change the number of elements)
        if numel(idx)~=this.nElems ; error('The provided indices must be of size [nElems 1]') ; end
        if any(idx>numel(this.Types)) ; error('The indices must range in the number of available element types') ; end 
        this.Indices = [idx(:) this.NodeIdx] ;
    end
    
    function this = set.NodeIdx(this,idx)
    % Set element node indices (CAN change the number of elements)
        if size(idx,1)==this.nElems % nothing to worry about
            this.Indices = [this.TypeIdx idx(:,:)] ;
        else % the number of elements change, clear the element types indices (that will be automatically detected)
            this.Indices = [zeros(size(idx,1),1,'uint32') idx(:,:)] ;
        end
    end
    
    function this = set.Types(this,types)
    % Set the list of elements
    % Change the type index list
        types = types(:) ;
        if this.nElems==0 % Do nothing, no elements in the list
            this.Types = types(:) ;
            return ;
        elseif numel(types)==this.nElems % One type by element, already ok
            typeIdx = 1:numel(types) ;
        else % We have to find which type fit with the given element node indices
            % Try to find the old types in the new types
                [~,indOldInNewTypes] = ismember(this.Types,types) ;
                valid = this.TypeIdx>0 ;
                valid(valid) = valid(valid) & indOldInNewTypes(this.TypeIdx(valid))>0 ;
            % Assign type indices related to found old types
                typeIdx = zeros(this.nElems,1,'uint32') ;
                typeIdx(valid) = indOldInNewTypes(this.TypeIdx(valid)) ;
            % Auto assign non-found types
                if any(~valid)
                    typeIdx(~valid) = this.assignElementTypeIdx(this.NodeIdx(~valid,:),types) ;
                end
        end
        this.Types = types ;
        this.TypeIdx = typeIdx ;
    end
    
    function this = set.Indices(this,indices)
    % Set the Indices property
        if isempty(indices) ; this.Indices = [] ; return ; end 
    % Check if element type indices are valid
        % In the right range
            valid = indices(:,1)>0 & indices(:,1)<=numel(this.Types) ;
        % The number of node indices is correct
            nNodes = reshape([this.Types.nNodes],[],1) ;
            valid(valid) = valid(valid) & sum(indices(valid,2:end)>0,2)==nNodes(indices(valid,1)) ;
    % Auto assign type indices if needed
        if any(~valid) 
            indices(~valid,1) = this.assignElementTypeIdx(indices(~valid,2:end)) ;
        end
    % Set
        this.Indices = indices ;
    end
    
    function [vars,ie] = getListOf(this,var,ie)
    % Dispatch elements properties and functions
    % "var" is a string or function handle
    % corresponding to a common element function or property
        if nargin<2 ; error('What do you want in the list ?') ; end
    % Extract a subtable if needed
        if nargin>2 ; this = this.subpart(ie) ; 
        else
            ie = 1:this.nElems ;
        end 
    % Apply to the element type list
        if ischar(var)
            if ismethod(this.Types,var)
                vars = this.Types.(var) ;
                if ~iscell(vars) ; vars = num2cell(vars,2:ndims(vars)) ; end
            else
                vars = {this.Types.(var)}  ;
            end
        elseif isa(var,'function_handle')
            vars = arrayfun(var,this.Types,'UniformOutput','false')  ;
        else ; error('Unsupported input') ; 
        end
    % Homogeneize if needed
        nd = cellfun(@numel,vars) ;
        vars = pkg.data.matchsize(2:max(nd),vars) ;
        n1 = cellfun(@size,vars,repmat({1},size(vars))) ;
    % Replicate
        vars = [{[]} reshape(vars,1,[])] ;
        vars = vars(this.TypeIdx+1) ;
    % Concatenate
        n1 = [0 ; n1(:)] ;
        ie = repelem(ie(:),n1(this.TypeIdx+1),1) ;
        vars = cat(1,vars{:}) ;
%         [vars,ii] = pkg.data.padcat(1,vars{:}) ;
%         ie = ie(ii) ;
    end
end


%% ELEMENT TABLE CLEANUP
methods
    function this = clean(this)
    % Clean the table
        if isempty(this) ; return ; end
        if isempty(this.Types) ; this = [] ; return ; end
    % Cull elements with invalid Type
        noType = ~(this.TypeIdx>0 & this.TypeIdx<=this.nTypes) ;
        if any(noType)
            this.Indices = this.Indices(~noType,:) ;
        end
    % Cull fully invalid elements
        fullZeros = ~any(this.NodeIdx,2) ;
        if any(fullZeros)
            this.Indices = this.Indices(~fullZeros,:) ;
        end
    % Cull elements with duplicated nodes
        dupNode = any(diff(double(sort(this.indicesWithNaNs,2)),1,2)==0,2) ;
        if any(dupNode)
            this.Indices = this.Indices(~dupNode,:) ;
        end
    % Cull elements with incorrect number of nodes
        nNodesInType = reshape([this.Types.nNodes],[],1) ;
        wrongNumber = this.nNodes~=nNodesInType(this.TypeIdx) ;
        if any(wrongNumber)
            this.Indices = this.Indices(~wrongNumber,:) ;
        end
    % Make the element list unique
        this = unique(this) ;
    end
end


%% CONNECTIVITY RETRIEVING
methods
    
    function M = sparse(this,values)
    % Return a sparse matrix containing values of different kinds
        if nargin<2 ; values = 'logical' ; end
        M = pkg.data.list2sparse(this.NodeIdx,values) ;
    end
    
    function M = contains(this,features)
    % Return a sparse matrix of logical values
    % M(i,j) = true if this(j,:) contains features(i,:)
        if features.nElems==0 || this.nElems==0
            M = sparse(features.nElems,this.nElems) ; 
            return ;
        end
    % The test uses the number of shared node indices
        nNodes = [0 features.Types.nNodes] ;
        nNodes = nNodes(features.TypeIdx+1) ;
    % Shared number of nodes between tables (sparse matrix)
        Mt = this.sparse('logical') ;
        Mf = features.sparse('logical') ;
        if size(Mt,1)~=size(Mf,1)
            sz = min(size(Mt,1),size(Mf,1)) ;
            M = Mf(1:sz,:)'*Mt(1:sz,:) ;
        else
            M = Mf'*Mt ;
        end
    % Nonzeros
        [ii,jj,vv] = find(M) ;
    % Contains
        bb = vv(:)==reshape(nNodes(ii),[],1) ;
        M = sparse(ii(bb),jj(bb),true,features.nElems,this.nElems) ;
    end
    
    function [table,elem2feat] = getTableOfUnique(this,features)
    % Return the table of unique element features (with duplicates)
    % featIdx returns the feature indices associated to each element in the
    % list
        [table,featElmt,featIdxForElmt] = getTableOf(this,features) ;
        [table,ie] = unique(table) ;
        if nargout<2 ; return ; end
        elem2feat = sparse(ie(:),featElmt(:),featIdxForElmt(:),table.nElems,this.nElems) ;
    end
    
    function [table,featElmt,featIdxInElmt] = getTableOf(this,features)
    % Return the complete table of element features (with duplicates)
    % Features are 'Faces' or 'Edges'
        if nargin<2 || ~ismember(features,{'Faces','Edges'})
            error('Wrong feature argument') ;
        end
    % Init feature table with types
        table = pkg.geometry.mesh.elements.ElementTable([this.Types.(features)]) ;
    % Feature hierarchy
        % Features counts
            nFeatInElemType = [0 this.Types.(['n' features])] ; % (includes typeIdx=0)
            nFeatTypesInPrevElems = cumsum(nFeatInElemType) ;
        % Element associated to each feature
            featElmt = repelem(1:this.nElems,nFeatInElemType(this.TypeIdx+1)) ;
        % Feature index in the parent element
            featIdxInElmt = arrayfun(@colon,nFeatInElemType*0 + 1,nFeatInElemType,'UniformOutput',false) ;
            featIdxInElmt = [featIdxInElmt{this.TypeIdx+1}] ;
        % Feature index in the current table
            featIdxInTable = arrayfun(@colon,nFeatTypesInPrevElems - nFeatInElemType + 1,nFeatTypesInPrevElems,'UniformOutput',false) ;
            featIdxInTable = [featIdxInTable{this.TypeIdx+1}] ;
    % "Local" node indices (in the ELEMET list of nodes & types)
        nodeIdx = table.NodeIdx(featIdxInTable,:) ;
    % "Global" node indices (in the MESH list of nodes)
        % Linear indexation
            valid = nodeIdx>0 ; % Valid feature indices
        	featElmtRep = repmat(featElmt(:),1,table.nMaxNodesByElem) ;
            iii = sub2ind(size(this.Indices),featElmtRep(valid),nodeIdx(valid)) ;
        % Global indices (MESH node number)
            nodeIdx(valid) = this.NodeIdx(iii) ;
    % "Global" type indices. Local typeIdx are forced to be valid (as defined in elements classes)
        typeIdx = table.TypeIdx(featIdxInElmt+nFeatTypesInPrevElems(this.TypeIdx(featElmt))) ;
    % Set the table
        table.Indices = [typeIdx(:) nodeIdx] ;
    end
end


%% DATA PROJECTION ON INDICES
methods
    function indices = indicesWithNaNs(this)
    % Return node indices with NaNs in place of zero's
        indices = double(this.NodeIdx) ;
        indices(indices<=0) = NaN ;
    end
    
    function VAL = dataAtIndices(this,DATA)
    % Return values VAL corresponding to the data DATA queried at this.NodeIdx
    % See pkg.data.dataAtIndices
        VAL = pkg.data.dataAtIndices(DATA,this.NodeIdx) ;
    end
    
    function VAL = meanDataAtIndices(this,DATA,DIM)
    % Return values VAL corresponding to the mean along the dimension DIM 
    % of data DATA queried at this.NodeIdx
    % See pkg.data.meanDataAtIndices
        if nargin<3
            VAL = pkg.data.meanDataAtIndices(DATA,this.NodeIdx) ;
        else
            VAL = pkg.data.meanDataAtIndices(DATA,this.NodeIdx,DIM) ;
        end
    end
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
    
    function [indices,iT] = catIndices(tables)
    % Concatenate indices of multiple tables
    % Retrieve indices
        IND = {tables.Indices} ;
    % Retrieve indices size
        [nElems,nMax] = cellfun(@size,IND) ;
        [uNMax,~,ia] = unique(nMax) ;
    % Corresponding input table number
        iT = repelem(1:numel(tables),nElems)' ;
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
    
    function [this,ia] = unique(this)
    % Return a table of UNIQUE elements
    % Make the list of element types unique
        [types,~,uTypeIdx] = unique(this.Types) ;
        validType = this.TypeIdx>0 ;
        this.TypeIdx(validType) = uTypeIdx(this.TypeIdx(validType)) ;
        this.Types = types ;
    % Cull duplicated element (same type AND node indices list)
        sortedIndices = [this.TypeIdx sort(this.NodeIdx,2)] ;
        [~,ie,ia] = unique(sortedIndices,'rows') ;
        if numel(ie)~=numel(ia)
            this.Indices = this.Indices(ie,:) ;
        end
    end
    
    function IDX = globalize(this,localIdx,elmtIdx)
    % Return global node indices (mesh level) corresponding to element
    % local node indices (reference element level)
    % localIdx [nElmtIdx nColumns] node indices in the reference element
    % elmtIdx (opt.) [nElmt 1] used to apply the method to a custom list of
    % elements (with some possibly repeated and/or removed)
        nElems = size(localIdx,1) ;
        if nargin<3 ; elmtIdx = 1:nElems ; end
        IDX = zeros(size(localIdx),'uint32') ;
    % Full global indices
        globIdx = this.NodeIdx(elmtIdx,:) ;
    % Find valid local indices
        [el,col,locIdx] = find(uint32(localIdx)) ;
    % Linear indexing
        IDX(sub2ind(size(IDX),el,col)) = globIdx(sub2ind(size(globIdx),el,locIdx)) ;
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
%     % Sub-table extraction with out = this(ind,:)
%         if strcmp(s(1).type,'()') % && numel(s(1).subs)==1
%             disp('TABLE SUBSREF') ;
%             if this.nElems==0 % No elements to extract
%                 varargout{1} = pkg.geometry.mesh.elements.ElementTable() ;
%             else % There is elements to extract 
%                 % Force extraction of all indices
%                     s(1).subs = {s(1).subs{1},':'} ; 
%                 % Extract the sub-list in indices
%                     indices = subsref(this.Indices,s(1)) ; 
%                 % List of used types
%                     [typeIdx,~,it] = unique(indices(:,1)) ;
%                     types = this.Types(typeIdx(typeIdx~=0)) ;
%                     indices(:,1) = it(:) ;
%                 % Create the sub-table
%                     varargout{1} = pkg.geometry.mesh.elements.ElementTable('Types',types,'Indices',indices) ;
%             end
%             if numel(s)>1 % If the subsref is deeper..
%                 [varargout{1:nargout}] = subsref(varargout{1},s(2:end)) ;
%             end
%         else
%             [varargout{1:nargout}] = builtin('subsref',this,s);
%         end
%     end

%     function this = subsasgn(this,s,data)
%     % Sub-assignment
%     % Allow to assign a sub-part of ONE table with element indices
%         if numel(this)>1 
%             [varargout{1:nargout}] = builtin('subsasgn',this,s) ; 
%             return ;
%         end
%     % Sub-table assignment with this(ind,:) = data
%         if strcmp(s(1).type,'()') % && numel(s(1).subs)==1
%             disp('TABLE SUBSASGN') ;
%             if numel(s)==1 % A sub-part of the table has to be set
%                 subTable = subsref(this,s(1)) ;
%             elseif numel(s)==2 % Some properties of a sub-part of the table has to be set
%             else % Cannot
%             end
%             % Force extraction of all indices
%                 s(1).subs = {s(1).subs{1},':'} ; 
%             % Extract the sub-list in indices
%                 indices = subsref(this.Indices,s(1)) ; 
%             % List of used types
%                 [typeIdx,~,it] = unique(indices(:,1)) ;
%                 types = this.Types(typeIdx(typeIdx~=0)) ;
%                 indices(:,1) = it(:) ;
%             % Create the sub-table
%                 varargout{1} = pkg.geometry.mesh.elements.ElementTable('Types',types,'Indices',indices) ;
%             if numel(s)>1 % If the subsref is deeper..
%                 [varargout{1:nargout}] = subsref(varargout{1},s(2:end)) ;
%             end
%         else % All other cases
%             [varargout{1:nargout}] = builtin('subsref',this,s);
%         end
%     end
    
%     function ii = end(this,k,~)
%     % Return the last indice in the table
%         switch k
%             case 1 ; ii = this.nElems ;
%             case 2 ; ii = size(this.Indices) ;
%             otherwise ; ii = 1 ;
%         end
%     end

    function table = subpart(this,elmtIdx)
    % Return a new table containing only a sub-part of the elements
    % Extract the sub-indices
        indices = this.Indices(elmtIdx,:) ;
        typeIdx = this.TypeIdx(elmtIdx) ;
    % Keep only the unique used types
        validType = typeIdx>0 ;
        [keepTypes,~,typeIdx(validType)] = unique(typeIdx(validType)) ;
        indices(:,1) = typeIdx ;
    % Assign
        table = pkg.geometry.mesh.elements.ElementTable('Types',this.Types(keepTypes),'Indices',indices) ;
    end
        
    
    function table = cat(~,varargin)
    % Concatenate element tables to an unique table
    % Convert to tables
        tables = cellfun(@pkg.geometry.mesh.elements.ElementTable,varargin) ;
    % Concatenate Element Types
        types = cat(1,tables.Types) ;
    % Concatenate Indices
        [indices,ta] = catIndices(tables) ;
        if isempty(indices)
            table = pkg.geometry.mesh.elements.ElementTable('Types',types) ; 
            return ; 
        end
    % Change Element Type Indices 
        nTypesInTable = [tables.nTypes] ;
        nTypesInPreviousTables = [0 cumsum(nTypesInTable)] ;
        typeIdx = indices(:,1) + uint32(reshape(nTypesInPreviousTables(ta),[],1)) ;
        typeIdx(indices(:,1)==0) = 0 ;
        indices(:,1) = typeIdx ;
    % Create the table
        table = pkg.geometry.mesh.elements.ElementTable() ; 
        table.Types = types ;
        table.Indices = indices ; 
    end
    function table = vertcat(varargin) ; table = cat(1,varargin{:}) ; end
    function table = horzcat(varargin) ; table = cat(1,varargin{:}) ; end
    
    function table = repmat(this,varargin)
    % Replicate vertically the element table
        if nargin==2 ; nRep = prod(varargin{1}) ; end
        if nargin>2 ; nRep = prod(cat(2,varargin{:})) ; end
        table = pkg.geometry.mesh.elements.ElementTable('Types',this.Types,'Indices',repmat(this.Indices,[nRep 1])) ;
    end
end

%% CONVERSION TO SIMPLEX MESH
methods
    function table = simplex(this)
    % Return a table representing the same mesh but with simplex elements
        table = [] ;
        for tt = 1:this.nTypes
            elmtTable = this.Types(tt).simplices ;
            nodeIdx = elmtTable.dataAtIndices(this.NodeIdx') ;
            nodeIdx = reshape(permute(nodeIdx,[1 3 2]),[],size(nodeIdx,2)) ;
            indices = padarray(nodeIdx,[0 1],1,'pre') ;
            newElems = pkg.geometry.mesh.elements.ElementTable('Types',elmtTable.Types,'Indices',indices) ;
            table = [table newElems] ;
        end
    end
end


%% ELEMENT IDENTIFICATION
methods
    function bool = isP1Quad(this)
    % Return true for each LagrangeElement('quad',1) in the table
        bool = this.Types.isP1Quad ;
        bool = ismember(this.TypeIdx,find(bool)) ;
    end
end


%% NODE LOCALIZATION ETC...
methods    
end

end
