classdef (Abstract) AbstractElement < handle & matlab.mixin.Heterogeneous
%ABSTRACTELEMENT Base class for mesh elements


%% (ABSTRACT) ELEMENT DEFINITION 
% ELEMENT GEOMETRY
properties (Abstract , SetAccess = protected)
    % Node local Coordinates [nNodes nDims]
    NodeLocalCoordinates
    % The list of faces (ElementTable) [nFaces nMaxNodesByFace]
    Faces pkg.geometry.mesh.elements.ElementTable
    % The list of edges (ElementTable) [nEdges 2] 
    Edges pkg.geometry.mesh.elements.ElementTable 
end
% SHAPE FUNCTIONS
methods (Abstract)
    % Evaluate the shape functions at local coordinates E
    % So that f(E) = N*f_nodes
    % E = [nE nDims] , N = [nE nNodes]
    N = evalAt(this,E)
end
% INTEGRATION QUADRATURE
properties (Abstract)
    GaussIntegrationPoints % [nGaussIntPts nDims]
    GaussIntegrationWeights % [nGaussIntPts 1]
end


%% GEOMETRIC INFORMATIONS
methods (Sealed)
    % Number of ...
    function val = nNodes(this) ; [val,~] = cellfun(@size,{this.NodeLocalCoordinates}) ; end
    function val = nFaces(this) ; val = cellfun(@nElems,{this.Faces}) ; end
    function val = nMaxNodesByFace(this) ; val = cellfun(@nMaxNodesByElem,{this.Faces}) ; end
    function val = nEdges(this) ; val = cellfun(@nElems,{this.Edges}) ; end
    function val = nMaxNodesByEdge(this) ; val = cellfun(@nMaxNodesByElem,{this.Edges}) ; end
    % Parameter space
    function val = nDims(this) ; [~,val] = cellfun(@size,{this.NodeLocalCoordinates}) ; end
    function val = localCoordinatesDomain(this) 
        val = [min(this.NodeLocalCoordinates,[],1) ; max(this.NodeLocalCoordinates,[],1)] ;
    end
    function val = localCoordinatesIn3D(this,e) 
        if nargin<2 ; e = this.NodeLocalCoordinates ; end
        val = [e zeros(size(e,1),3-size(e,2))] ;
    end
    % Circumcenter
    function C = centroid(this) ; C = permute(mean(cat(3,this.NodeLocalCoordinates),1),[3 2 1]) ; end
end

%% SHAPE FUNCTIONS DERIVATIVES
methods

    function delta = defaultDelta(this)
    % Default step used for finite diff gradient evaluation
        delta = range(this.localCoordinatesDomain,1)*1e-6 ;
    end

    function DER = evalDerivativeAt(this,E,ORD,delta)
    % Evaluate the partial derivatives of order ORD = [o1 o2 ...] (:[1 nDims]) of 
    % shape functions at given local coordinates E:[nRows nDims]
    % Evaluation performed by recursive finite differences. 
    % When possible, it is preferred to override this method to provide
    % analytical derivatives !
    % exemple: df(E)/(dx²dy) = DER(E,[2 1 0]) ;
        if nargin<3 ; ORD = [1 zeros(1,this.nDims-1)] ; end
        if nargin<4 ; delta = this.defaultDelta ; end
        if numel(ORD)~=this.nDims ; error('Wrong derivation order argument (must be [1 nDims])') ; end
        % Evaluation by a centered 2nd-order finite diff. scheme
            if all(ORD==0) % No Derivative
                DER = this.evalAt(E) ;
            else % Recursive call
                DIM = find(ORD~=0,1,'first') ;
                newORD = ORD-sparse(1,DIM,1,1,this.nDims) ;
                de = (ORD-newORD).*delta ;
                DER = ( this.evalDerivativeAt(E+de,newORD,delta) ...
                        - this.evalDerivativeAt(E-de,newORD,delta) ...
                        ) / ( 2 * delta(DIM) ) ;
            end
    end

    function J = evalJacobianAt(this,E,delta)
    % Return the shape function gradient evaluated at local coordinates E
    % E = [nE nDims]
    % delta = optionnal, for finite diff derivatives (see above)
    % J = [nE nNodes nDims]
    % So that dN_dei = J(:,:,i)*Ne
        if nargin<3 ; delta = this.defaultDelta ; end
        J = NaN(size(E,1),this.nNodes,this.nDims) ;
        for dim = 1:this.nDims
            ORD = full(sparse(1,dim,1,1,this.nDims)) ;
            J(:,:,dim) = this.evalDerivativeAt(E,ORD,delta) ;
        end
    end

    function [H,ORD,id] = evalHessianAt(this,E,delta)
    % Return the shape function second derivative evaluated at local coordinates E
    % E = [nE nDims]
    % delta = optionnal, for finite diff derivatives (see above)
    % G = [nE nNodes nDims nDims]
    % So that d²N_(dei.dej) = H(:,:,i,j)*Ne
    % If only H is queried, then H is the hessian matrix (with
    % symmetric values duplicated)
    % Otherwise, ORD is the orders of derivations corresponding to each
    % slice of H and id contains the index of the slice of H
    % corresponding to the right derivative: Hessian = H(:,:,id) ;
        if nargin<3 ; delta = this.defaultDelta ; end
    % Build the derivation order vectors
        ORD = reshape(eye(this.nDims),[1 this.nDims this.nDims]) ;
        ORD = ORD + permute(ORD,[2 1 3]) ;
        ORD = reshape(ORD,[],this.nDims) ;
        [ORD,~,id] = unique(ORD,'rows') ; % prevent computation of the same derivatives
        id = reshape(id,[this.nDims this.nDims]) ;
        nDeriv = size(ORD,1) ;
    % Intialize the Hessian
        H = NaN(size(E,1),this.nNodes,nDeriv) ;
    % Compute
        for dd = 1:nDeriv
            H(:,:,dd) = this.evalDerivativeAt(E,ORD(dd,:),delta) ;
        end
    % If the hessian matrix is asked for
        if nargout==1 ; H = reshape(H(:,:,id),size(E,1),this.nNodes,this.nDims,this.nDims) ; end
    end
end


%% GEOMETRY (TESTS and EVALUATION)
methods
    function tol = defaultTolerance(this)
    % Default tolerance for geometrical tests
        tol = norm(range(this.NodeLocalCoordinates,1))*1e-9 ;
    end
    
    function table = elementTable(this)
    % Return an ElementTable corresponding to the element alone
        table = pkg.geometry.mesh.elements.ElementTable(...
                    'Types',this ...
                    ,'Indices',[1 1:this.nNodes]) ;
    end
    
    function m = mesh(this,X)
    % Return a mesh corresponding to the element
        if nargin<2 ; X = this.NodeLocalCoordinates ; end
        table = elementTable(this) ;
        m = pkg.geometry.mesh.Mesh(...
                'Nodes',X ...
                ,'Elems',table) ;
    end
    
    function [in,on] = isInside(this,E,tol)
    % Return true if the local coordinates E are in/on the reference element
    % Uses the fact that the reference element is (SOULD BE) convex !
    % AND the fact that the edges and faces are anti-cockwise oriented !
    % AND reference faces/edges are supposed flat/straight
    % Should be overriden by inherited elements to provide an analytical
    % formulation
    % E = [nE nDims]
    % (in,on) = [nE 1] logical
        if nargin<3 ; tol = this.defaultTolerance ; end
    % Test the element bounding box first
        bbox = localCoordinatesDomain(this) ;
        in = all(E>=bbox(1,:)-tol & E<=bbox(2,:)+tol,2) ;
        on = in & all(E<=bbox(1,:)+tol | E>=bbox(2,:)-tol,2) ;
        if ~any(in) ; return ; end
    % 1D elements: just inside the domain==bbox
        if this.nDims==1 ; return ;
        else
    % 2D elements: all scal(vect(edgecenters->E),edgenormals)<0
        if this.nDims==2 ; features = this.Edges ;
    % 3D CONVEX elements: all scal(vect(facecenters->E),facenormals)<0
        else ; features = this.Faces ; 
    % Process 
        end
            nodes = features.dataAtIndices(this.NodeLocalCoordinates) ; % [nFeatures nMaxNodesByFeature nDims]
            midNode = permute(mean(nodes,2,'omitnan'),[2 3 1]) ; % [1 nDims nFeatures]
            if this.nDims==2 % 2D elements -> edge normals
                tangent = permute(diff(nodes(:,1:2,:),1,2),[2 3 1]) ; % Tangent vector [1 nDims nFeatures]
                normals = [tangent(:,2,:) -tangent(:,1,:)] ; % 90° clockwise rotation
            else % 3D elements -> Face normals
                T = zeros(this.nDims,2,features.nElems) ;
                for ff = 1:features.nElems
                    face = features.Types(features.TypeIdx(ff)) ;
                    dN_de = face.evalJacobianAt(face.centroid) ; % [1 face.nNodes face.nDims(==2)]
                    dx_de = sum(dN_de.*permute(nodes(ff,1:face.nNodes,:),[3 2 1]),2) ; % [this.nDims 1 face.nDims(==2)] 
                    T(:,:,ff) = squeeze(dx_de) ; % [nDims 2 1]
                end
                normals = [T(2,1,:).*T(3,2,:)-T(3,1,:).*T(2,2,:) T(3,1,:).*T(1,2,:)-T(1,1,:).*T(3,2,:) T(1,1,:).*T(2,2,:)-T(2,1,:).*T(1,2,:)] ; % [1 nDims nFeatures]
            end
            vect = E(in,:)-midNode ; % vect(edgecenters->E) [nE nDims nFeatures]
            scal = sum(vect.*normals,2) ; % scalar product [nE 1 nFeatures]
            in0 = in ;
            in(in0) = all(scal<=tol,3) ;
            on(in0) = in(in0) & any(abs(scal)<=tol,3) ;
        end
    end
    
    function F = evalFrameAt(this,E)
    % Return the element 3D frame matrix
    % so that frame(E) = F(E)*X
    end
    
    function [elems,idx] = slice(this,nodeBool)
    % Slice the element given a signed logical value of a levelset on
    % each node: -1 (inside), 0 (on) or 1 (outside)
    % input: this: element; nodeBool [nElems this.nNodes]
    % output:
    %   - an element table ELEMS containing the sliced elements. 
    %       /!\ the node indices in the table are complex uint32 ! :
    %       - real indices denote nodes of the reference element
    %       - imaginary indices denote edges of the reference element
    %   - a list IDX of size [ELEMS.nElems 1] where IDX(i) contains the
    %   index of the input element
        elems = pkg.geometry.mesh.elements.ElementTable ;
        idx = [] ;
    end
        
    function S = sliceCases(this)
    % Return all the cases that have to be tested while slicing this
    % element
    % S is an array of structures with fields:
    %   - Test: signed boolean configuration [1 this.nNodes]
    %           (-1 inside, 0 on, 1 outside, NaN any cases)
    %   - IN: elements inside (pkg.geometry.mesh.elements.ElementTable)
    %   - ON: elements on the slice (pkg.geometry.mesh.elements.ElementTable)
    %   - OUT: elements outside (pkg.geometry.mesh.elements.ElementTable)
    % /!\ The parent function will already test opposite signs, do not bother with this !
    % in the ElementTables, indices>this.nNodes denote edges indices
    % (used when the slice cross an edge)
        import pkg.geometry.mesh.elements.ElementTable
        import pkg.geometry.mesh.elements.base.*
        S = repmat(struct('Test',[],'IN',[],'ON',[],'OUT',[]),[0 0]) ; % empty structure
    % By default, this will return all test possibilities...
        % Create all possible combinations...
            allTests = repmat({[-1,1]},[this.nNodes 1]) ;
            [allTests{:}] = ndgrid(allTests{:}) ;
            allTests = cat(this.nNodes+1,allTests{:}) ;
            allTests = reshape(allTests,[],this.nNodes) ;
        % Sort by absolute value
            [~,ind] = sort(sum(abs(allTests),2)) ;
            allTests = allTests(ind,:) ;
        % Inject in the structure
            allTests = num2cell(allTests,2) ;
            [S(1:numel(allTests)).Test] = deal(allTests{:}) ;
    end
end


%% CONVERT TO SIMPLICES
methods
    function table = simplices(this)
    % Return a table of simplices representing the same element
        if this.nDims<=1 ; table = elementTable(this) ; return ; end
        nodeIdx = delaunay(this.NodeLocalCoordinates) ;
        indices = padarray(nodeIdx,[0 1],1,'pre') ;
        switch this.nDims
            case 2
                type = pkg.geometry.mesh.elements.base.Triangle ;
            case 3
                type = pkg.geometry.mesh.elements.base.Tetrahedron ;
        end
        table = pkg.geometry.mesh.elements.ElementTable('Types',type,'Indices',indices) ;
    end
    
    function mesh = simplexMesh(this)
    % Return a mesh of simplices representing the same element
        mesh = pkg.geometry.mesh.Mesh('Nodes',this.NodeLocalCoordinates,'Elems',this.simplices) ;
    end
end


%% PLOT FUNCTIONS
methods
    function H = plot(this,toPlot,varargin)
    % Handles plot functions
        if nargin<2 ; toPlot = 'ReferenceElement' ; end
        if nargin<3 ; varargin = {} ; end
        switch upper(toPlot)
            case 'REFERENCEELEMENT'
                H = plotReferenceElement(this,varargin{:}) ;
            case 'SHAPEFUNCTIONS'
                H = plotShapeFunctions(this,varargin{:}) ;
            case 'SLICECASES'
                H = plotSliceCases(this,varargin{:}) ;
        end
    end

    function H = plotReferenceElement(this,ax)
    % Plot the reference element in a given axes
        if nargin<2 ; ax = gca ; end
        h = gobjects(0) ;
        mesh = this.mesh(this.localCoordinatesIn3D) ;
        H = plot(mesh) ;
        H.Faces.FaceColor = 'none' ;
        H.VisibleNodes = 'all' ;
        H.HighlightEndNodes = false ;
        H.ShowLabels = {'Nodes','Edges','Faces'} ;
        if this.nDims>2 ; set(ax,'view',[30 30],'Toolbar',[],'Interactions',[rotateInteraction]) ; end
    end

    function H = plotShapeFunctions(this,fig)
    % Plot the shape functions in the current figure
        if nargin<2 ; fig = gcf ; end
        subdiv = 25*ones(1,this.nDims) ; % Number of Subdivisions
        switch this.nDims % Displayed derivation orders
            case 1
                ORD = [0;1;2] ;
            case 2
                ORD = [0 0 ; 1 0 ; 0 1 ; 2 0 ; 0 2 ; 1 1] ;
            case 3
                ORD = [0 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1] ;
        end
        % Local coordinates
            bbox = this.localCoordinatesDomain ;
            E = arrayfun(@(dim)linspace(bbox(1,dim),bbox(2,dim),subdiv(dim)),1:this.nDims,'UniformOutput',false) ;
            if numel(E)>1 ; [E{:}] = ndgrid(E{:}) ; end
            EE = cat(this.nDims+1,E{:}) ;
        % Delete points outside the element
            e = reshape(EE,[],this.nDims) ;
            e(~isInside(this,e)) = NaN ;
            EE = reshape(e,size(EE)) ;
        % Figure
            clf(fig) ;
            H = gobjects(0) ;
            ax = gobjects(0) ;
            lin = gobjects(0) ;
            srf = gobjects(0) ;
            scat = gobjects(0) ;
            ttl = gobjects(0) ;
         % For each derivation order
            nOrd = size(ORD,1) ;
            for oo = 1:nOrd
                % Evaluate the shape function derivatives
                    N = this.evalDerivativeAt(e,ORD(oo,:)) ;
                    NN = reshape(N,[subdiv this.nNodes]) ;
                % For each shape function
                    for nn = 1:this.nNodes
                        % Init axes
                            ax(end+1) = axes ;
                            ax(end).OuterPosition = [(nn-1)/this.nNodes 1-oo/nOrd 1/this.nNodes 1/nOrd] ;
                        % Plot the reference element
                            H = this.plot('ReferenceElement',ax(end)) ;
                            H.ShowLabels = {} ;
                        % Remove the labels
                            H.ShowLabels = {} ;
                        % Highlight the current node only
                            H.VisibleNodes = 'all' ;
                            H.HighlightEndNodes = false ;
                            H.Nodes.Faces = nn ;
                            H.Nodes.MarkerSize = 30 ;
                        % Plot the shape function/derivative
                            switch this.nDims
                                case 1
                                    lin(end+1) = plot(ax(end),EE(:,:,1),NN(:,nn)) ;
                                case 2
                                    srf(end+1) = surf(ax(end),EE(:,:,1),EE(:,:,2),NN(:,:,nn),NN(:,:,nn)) ;
                                case 3
                                    scat(end+1) = scatter3(e(:,1),e(:,2),e(:,3),10,N(:,nn)) ;
                                otherwise % ???
                            end
                        % Add a title
                            ttl(end+1) = title(ax(end),['N_{' num2str(nn) '}']) ;
                            if any(ORD(oo,:)) 
                                ttl(end).String = [ttl(end).String(1:end-1) ...
                                                    ',' sprintf('%i',repelem(1:this.nDims,1,ORD(oo,:))) ...
                                                    ttl(end).String(end)] ; 
                            end 
                    end
            end
        % Set Common properties
            set(ax,'FontSize',12,'TickLabelInterpreter','none') ;
            if this.nDims<2
                set(ax,'view',[0 90],'Toolbar',[],'Interactions',[]) ;
            else
                set(ax,'view',[30 30],'Toolbar',[],'Interactions',[rotateInteraction]) ;
                set(srf,'FaceColor','interp','EdgeColor','none')
                fig.UserData.hlink = linkprop(ax,'view') ;
            colormap(jet(11)) ;
            end
            set(ttl,'Interpreter','tex','units','normalized','position',[0.5 0.5 0]) ;
    end
    
    function ax = plotSliceCases(this,fig,cases)
    % Plot the slice cases implemented for the element
        if nargin<2 ; fig = gcf ; end
        S = this.sliceCases ;
        if nargin<3 ; cases = 1:min(numel(S),30) ; end
        if max(cases)>numel(S) ; error('Invalid cases queried.') ; end
        nCases = numel(cases) ;
        nAx = ceil(sqrt(nCases)) ; nAy = ceil(nCases/nAx) ;
        ax = gobjects(0) ;
        E = this.localCoordinatesIn3D ;
        E = [E ; this.Edges.meanDataAtIndices(E)] ; % add edge centroids
        clrmp = padarray(linspecer(3),10,'replicate') ; jet(101) ;
        refMesh = this.mesh ;
        for cc = cases(:)'
        % New axes
            nn = numel(ax)+1 ;
            yy = ceil(nn/nAx) ; xx = nn-(yy-1)*nAx ;
            ax(end+1) = axes(fig,'OuterPosition',[(xx-1)/nAx 1-yy/nAy 1/nAx 1/nAy]) ; 
        % Reference element
            refMesh.Nodes = E(1:this.nNodes,:) ;
            H = plot(refMesh) ;
        % Show the levelset values
            H.VisibleFaces = 'all' ;
            H.Faces.FaceColor = 'interp' ;
            H.Faces.FaceVertexCData = reshape(S(cc).Test(1:this.nNodes),[],1) ;
        % Show inside, onside, outside meshes
            parts = {'IN','ON','OUT'} ;
            shiftE = 1.1*[1 0 0 ; 1 -1 0 ; 0 -1 0] ;
            partColors = clrmp([1 ceil(end/2) end],:) ;
            faceAlpha = 0.5 ;
            for pp = 1:numel(parts)
                refMesh.Nodes = E(1:this.nNodes,:) + shiftE(pp,:) ;
                h = plot(refMesh,'VisibleFaces','none','HighlightBoundaryEdges',false) ;
                    h.ShowLabels = {} ;
                    h.Faces.FaceColor = 'none' ;
                    h.HighlightBoundaryEdges = false ;
                partMesh = pkg.geometry.mesh.Mesh('Nodes',E + shiftE(pp,:),'Elems',S(cc).(parts{pp})) ;
                if partMesh.nElems>0
                    Hp = plot(partMesh) ;%,'VisibleNodes','all') ;
                        Hp.Faces.FaceColor = partColors(pp,:) ;
                        Hp.Faces.FaceAlpha = faceAlpha ;
                        Hp.Edges.EdgeColor = partColors(pp,:) ;
                        %Hp.Nodes.MarkerEdgeColor = partColors(pp,:) ;
                        Hp.BoundaryEdges.EdgeColor = partColors(pp,:) ;
                end
            end
        end
        set(ax,'xtick',[],'ytick',[],'ztick',[]) ;
        set(ax,'looseinset',[1 1 1 1]*0.05) ;
        set(ax,'colormap',clrmp) ;
        set(ax,'clim',[-1 1]) ;
        axis(ax,'off') ;
        axis(ax,'equal') ;
        fig.UserData.hlink = linkprop(ax,'view') ;
    end
end

%% ELEMENT IDENTIFICATION FUNCTION
methods (Sealed)
    function bool = isP1Quad(elems)
    % Return true for LagrangeElements('quad',1)
    % Useful to switch node numbering from trigo to P1 quad nodes
        elems = elems(:)' ;
        bool = [elems.nNodes] == 4 ; 
        if any(bool) ; bool = bool & arrayfun(@(e)isa(e,'pkg.geometry.mesh.elements.LagrangeElement'),elems) ; end
        if any(bool) ; bool(bool) = bool(bool) & [elems(bool).Order]==1 ; end
        if any(bool) ; bool(bool) = bool(bool) & strcmp({elems(bool).Geometry},'quad') ; end
    end
end

end

