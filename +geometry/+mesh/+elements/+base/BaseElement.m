classdef (Abstract) BaseElement < pkg.geometry.mesh.elements.AbstractElement
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
        
    function S = sliceCases(this)
    % Return all the cases that have to be tested while slicing this
    % element
    % S is an array of structures with fields:
    %   - Test: signed boolean configuration [1 this.nNodes]
    %           (-1 inside, 0 on, 1 outside, NaN any cases)
    %   - IN: elements inside (pkg.geometry.mesh.elements.ElementTable)
    %   - ON: elements on the slice (pkg.geometry.mesh.elements.ElementTable)
    %   - OUT: elements outside (pkg.geometry.mesh.elements.ElementTable)
    % in the ElementTables, NodeIdx>this.nNodes denote edges indices
    % (used when the slice cross an edge)
        import pkg.geometry.mesh.elements.ElementTable
        S = repmat(struct('Test',[],'IN',[],'ON',[],'OUT',[]),[0 0]) ; % empty structure
    % Element Types
        import pkg.geometry.mesh.elements.base.*
        node = Node ; bar = Bar ; tri = Triangle ; tet = Tetrahedron ;
    % By default, this will return all test possibilities... 
    % Create all possible combinations (more than 6000 tests for an hexahedron...)
        nodeTests = pkg.math.combinations(repmat({[-1,1]},[this.nNodes 1])) ;
    % Corresponding crossing edges
    % an edge cross if the sign of the lvlset is opposite on the two ends
        crossEdg = this.Edges.dataAtIndices(nodeTests') ;
        crossEdg = reshape(range(crossEdg,2)==2,this.nEdges,size(nodeTests,1)).' ;
        tests = [nodeTests crossEdg] ;
    % Inject tests in the structure
        testsCells = num2cell(tests,2) ;
        [S(1:numel(testsCells)).Test] = deal(testsCells{:}) ;
    % Trivial cases: the element is in/on/out -> all nodes have the same signe boolean)
        elem = ElementTable('Types',this,'Indices',[1 1:this.nNodes]) ;
        [S(all(nodeTests==-1,2)).IN] = deal(elem) ;
        [S(all(nodeTests==0,2)).ON] = deal(elem) ;
        [S(all(nodeTests==1,2)).OUT] = deal(elem) ;
    % Non-trivial cases: we create a mesh of simplices from nodes with the same status
        nonTrivial = find(range(nodeTests,2)>0) ; nNT = numel(nonTrivial) ;
        if nNT==0 ; return ; end
        % Logical indexing
            nodON = nodeTests(nonTrivial,:)==0 ;
            nodIN = nodeTests(nonTrivial,:)<0 ;
            nodOUT = nodeTests(nonTrivial,:)>0 ;
            edgCross = crossEdg(nonTrivial,:) ;
            idxAll = [true(nNT,this.nNodes) edgCross] ;
        % Simplex mesh
            if this.nDims==1 % 1D BAR elements
            % Just make a mesh from the two end nodes and the central point
                meshes = repmat({[1 3 ; 3 2]},[nNT 1]) ;
            else
            % Node coordinates for the delaunay mesh
                Ne = this.NodeLocalCoordinates ;
                Ne = [Ne ; this.Edges.meanDataAtIndices(Ne)] ; % edd "edge centroids"
            % Delaunay mesh corresponding to unique cases
                [~,ia,ic] = unique(edgCross,'rows') ;
                NE = cellfun(@(ii)Ne(ii,:),num2cell(idxAll(ia,:),2),'UniformOutput',false) ;
                meshes = cellfun(@delaunay,NE,'UniformOutput',false) ;
            % Spread on non-unique cases
                meshes = meshes(ic) ;
            end
    % Create the element tables
        for cc = 1:nNT
        % Retrieve global indices
            nodIdx = find(idxAll(cc,:)) ;
            mesh = nodIdx(meshes{cc}) ;
        % Get elems belonging to the different parts
        % Due to the nature of delaunay triangulation, an element cannot
        % contain both inside nodes and outside nodes
            elmtIN = any(ismember(mesh,find(nodIN(cc,:))),2) ;
            elmtOUT = any(ismember(mesh,find(nodOUT(cc,:))),2) ;
        % Fill
            switch this.nDims
                case 1 
                    S(nonTrivial(cc)).ON = ElementTable('Types',node,'Indices',[1 3]) ;
                    S(nonTrivial(cc)).IN = ElementTable('Types',bar,'Indices',[1 mesh(elmtIN,:)]) ;
                    S(nonTrivial(cc)).OUT = ElementTable('Types',bar,'Indices',[1 mesh(elmtOUT,:)]) ;
                case 2 
                % Build edges
                    edges = reshape(mesh(:,[1 2 2 3 3 1])',2,[])' ;
                    [edges,~,ec] = unique(sort(edges,2),'rows') ;
                    elem2edge = sparse(ec,repelem(1:size(mesh,1),3),true) ;
                % Get edges ON the levelset (all nodes are ON and it is
                % linked to at least one elmtIN or one elmtOUT)
                    edgON = all(ismember(edges,find([nodON(cc,:) edgCross(cc,:)])),2) ;
                    edgON = edgON(:) & any(elem2edge.*(elmtIN | elmtOUT)',2) ;
                % Fill the structure
                    S(nonTrivial(cc)).ON = ElementTable('Types',bar ...
                                    ,'Indices',padarray( edges(edgON,:) ,[0 1],1,'pre') ) ;
                    S(nonTrivial(cc)).IN = ElementTable('Types',tri ...
                                    ,'Indices',padarray( mesh(elmtIN,:) ,[0 1],1,'pre') ) ;
                    S(nonTrivial(cc)).OUT = ElementTable('Types',tri ...
                                    ,'Indices',padarray( mesh(elmtOUT,:) ,[0 1],1,'pre') ) ;
                case 3
                % Build faces
                    faces = reshape(mesh(:,[1 3 2 1 2 4 2 3 4 3 1 4])',3,[])' ;
                    [faces,~,fc] = unique(sort(faces,2),'rows') ;
                    elem2face = sparse(fc,repelem(1:size(mesh,1),4),true) ;
                % Get faces ON the levelset (all nodes are ON and it is
                % linked to at least one elmtIN or one elmtOUT)
                    faceON = all(ismember(faces,find([nodON(cc,:) edgCross(cc,:)])),2) ;
                    faceON = faceON(:) & any(elem2face.*(elmtIN | elmtOUT)',2) ;
                % Fill the structure
                    S(nonTrivial(cc)).ON = ElementTable('Types',tri ...
                                    ,'Indices',padarray( faces(faceON,:) ,[0 1],1,'pre') ) ;
                    S(nonTrivial(cc)).IN = ElementTable('Types',tet ...
                                    ,'Indices',padarray( mesh(elmtIN,:) ,[0 1],1,'pre') ) ;
                    S(nonTrivial(cc)).OUT = ElementTable('Types',tet ...
                                    ,'Indices',padarray( mesh(elmtOUT,:) ,[0 1],1,'pre') ) ;
            end
        end
    end
end

end

