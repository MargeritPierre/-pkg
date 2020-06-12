classdef (Abstract) BaseElement < pkg.mesh.elements.AbstractElement
%BASEELEMENT Abstract superclass for BASE mesh elements

%% LOGICAL TESTS & SORTING
methods (Sealed)
    function [uTypes,ia,ic] = unique(types,varargin)
    % Used to return a unique list of types
        classes = arrayfun(@class,types,'UniformOutput',false) ;
        [~,ia,ic] = unique(classes,varargin{:}) ;
        uTypes = types(ia) ;
    end
    
    function val = eq(a,b)
        classA = arrayfun(@class,a,'UniformOutput',false) ;
        classB = arrayfun(@class,b,'UniformOutput',false) ;
        val = reshape(strcmp(classA,classB),size(a)) ;
    end
end

%% ELEMENT CUTTING & SLICING
methods
    function S = cut(this,nodeBool) 
    % Cut an element with the given sign of a levelset at nodes
    % input: nodeBool [nCases this.nNodes] (-1) inside or (1) outside
    %   cases where the node is ON the levelset are not implemented
    % output: structure S with 6 fields:
    %   - elems(ON/IN/OUT): Table of elements containing the sliced element
    %       connectivities and types
    %       /!\ any elems.NodeIdx>this.nNodes denotes an edge index: 
    %       EdgIdx = NodeIdx-this.nNodes
    %   - idx(ON/IN/OUT): index of the case in nodeBool corresponding to
    %       each output element. 
    %       Used for node index globalization (see pkg.mesh.Mesh.cut)
        import pkg.mesh.elements.ElementTable
        import pkg.mesh.elements.base.*
        S = struct('elemsON',[],'idxON',[],'elemsIN',[],'idxIN',[],'elemsOUT',[],'idxOUT',[]) ;
        nCases = size(nodeBool,1) ;
        if nElmts==0 ; return ; end
    % Get edges crossing the level set (end nodes with two different sign)
        eBool = this.Edges.dataAtIndices(nodeBool') ;
        edgCross = reshape(eBool(:,1,:)~=eBool(:,2,:),[this.nEdges nCases]).' ;
    % 1D elements are straightforward
        if this.nDims==1
            S.idxON = find(edgCross(:,1)) ;
            S.elemsON = ElementTable('Types',Node,'Indices',[1 1])
            return ;
        end
    % Possible node coordinates
        Ne = this.NodeLocalCoordinates ;
        Ne = [Ne ; this.Edges.meanDataAtIndices(Ne)] ;
    % Process each case: create a delaunay triangulation then split it
    % for nDims>1, only an even number of edges can cross the levelset (
        for cc = 1:nCases
            
        end
    end
        
    function S = sliceCases(this)
    % Return all the cases that have to be tested while slicing this
    % element
    % S is an array of structures with fields:
    %   - Test: signed boolean configuration [1 this.nNodes]
    %           (-1 inside, 0 on, 1 outside, NaN any cases)
    %   - IN: elements inside (pkg.mesh.elements.ElementTable)
    %   - ON: elements on the slice (pkg.mesh.elements.ElementTable)
    %   - OUT: elements outside (pkg.mesh.elements.ElementTable)
    % in the ElementTables, NodeIdx>this.nNodes denote edges indices
    % (used when the slice cross an edge)
        import pkg.mesh.elements.ElementTable
        S = repmat(struct('Test',[],'IN',[],'ON',[],'OUT',[]),[0 0]) ; % empty structure
    % By default, this will return all test possibilities... 
    % (more than 6000 tests for an hexahedron...)
        % Create all possible combinations...
            allTests = repmat({[-1,1]},[this.nNodes 1]) ;
            [allTests{:}] = ndgrid(allTests{:}) ;
            allTests = cat(this.nNodes+1,allTests{:}) ;
            allTests = reshape(allTests,[],this.nNodes) ;
        % Sort by absolute value
            [~,ind] = sort(sum(abs(allTests),2)) ;
            allTests = allTests(ind,:) ;
        % Inject in the structure
            allTestsCells = num2cell(allTests,2) ;
            [S(1:numel(allTestsCells)).Test] = deal(allTestsCells{:}) ;
    % Corresponding crossing edges
    % an edge cross if the sign of the lvlset is opposite on the two ends
        crossEdg = this.Edges.dataAtIndices(allTests') ;
        crossEdg = reshape(range(crossEdg,2)==2,this.nEdges,numel(S)).' ;
    % Node coordinates for the delaunay mesh
        Ne = this.NodeLocalCoordinates ;
        Ne = [Ne ; this.Edges.meanDataAtIndices(Ne)] ;
    % Element Types
        import pkg.mesh.elements.base.*
        node = Node ; bar = Bar ; tri = Triangle ; tet = Tetrahedron ;
    % Create the element tables
        for tt = 1:numel(S)
            nodTest = S(tt).Test ;
            edgTest = find(crossEdg(tt,:)) + this.nNodes ;
        % Trivial cases: the full element is..
            if range(nodTest)==0 % all signs are the same
                elem = ElementTable('Types',this,'Indices',[1 1:this.nNodes]) ;
                switch nodTest(1)
                    case 0  ; S(tt).ON = elem ;
                    case -1 ; S(tt).IN = elem ;
                    otherwise ; S(tt).OUT = elem ;
                end
                continue ;
            end
        % Non-trivial cases: use delaunay triangulation...
        % Get the indices corresponding to each part
            idxON = [find(nodTest==0) edgTest] ;
            idxIN = [find(nodTest<=0) edgTest] ;
            idxOUT = [find(nodTest>=0) edgTest] ;
            idxAll = [1:this.nNodes edgTest] ;
        % Create a mesh of simplices from the usable node coordinates
            switch this.nDims
                case 1
                    allElems = [1 3 ; 3 2]' ;
                    allNodes = [1 2 3] ;
                case 2
                    allElems = idxAll(delaunay(Ne(idxAll,:))) ;
                    allEdges = reshape(allElems(:,[1 2 2 3 3 1])',2,[])' ;
                case 3
                    allElems = idxAll(delaunay(Ne(idxAll,:))) ;
                    allFaces = reshape(allElems(:,[1 3 2 1 2 4 2 3 4 3 1 4])',3,[])' ;
            end
        % Get elems belonging to the different parts
            elIN = allElems(all(ismember(allElems,idxIN),2),:) ;
            elOUT = allElems(all(ismember(allElems,idxOUT),2),:) ;
        % Fill
            switch this.nDims
                case 1 
                    S(tt).ON = ElementTable('Types',node,'Indices',[1 allNodes(idxON)]) ;
                    S(tt).IN = ElementTable('Types',bar,'Indices',[1 elIN]) ;
                    S(tt).OUT = ElementTable('Types',bar,'Indices',[1 elOUT]) ;
                case 2 
                    S(tt).ON = ElementTable('Types',bar ...
                                    ,'Indices',padarray( allEdges( all(ismember(allEdges,idxIN),2) ,:) ,[0 1],1,'pre') ) ;
                    S(tt).IN = ElementTable('Types',tri ...
                                    ,'Indices',padarray(elIN ,[0 1],1,'pre') ) ;
                    S(tt).OUT = ElementTable('Types',tri ...
                                    ,'Indices',padarray(elOUT ,[0 1],1,'pre') ) ;
                case 3
            end
        end
    end
end

end

