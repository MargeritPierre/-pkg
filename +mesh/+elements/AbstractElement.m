classdef (Abstract) AbstractElement < handle & matlab.mixin.Heterogeneous
%ABSTRACTELEMENT Base class for mesh elements


%% (ABSTRACT) ELEMENT DEFINITION 
% ELEMENT GEOMETRY
properties (Abstract , SetAccess = protected)
    % Node local Coordinates [nNodes nDims]
    NodeLocalCoordinates
    % The list of faces (ElementTable) [nFaces nMaxNodesByFace]
    Faces pkg.mesh.elements.ElementTable
    % The list of edges (ElementTable) [nEdges 2] 
    Edges pkg.mesh.elements.ElementTable 
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
    function C = circumcenter(this) ; C = permute(mean(cat(3,this.NodeLocalCoordinates),1),[3 2 1]) ; end
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
    % exemple: df(E)/(dx�dy) = DER(E,[2 1 0]) ;
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
    % G = [nE nNodes nDims]
    % So that dN_dei = G(:,:,i)*Ne
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
    % So that d�N_(dei.dej) = H(:,:,i,j)*Ne
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
                normals = [tangent(:,2,:) -tangent(:,1,:)] ; % 90� clockwise rotation
            else % 3D elements -> Face normals
                T = zeros(this.nDims,2,features.nElems) ;
                for ff = 1:features.nElems
                    face = features.Types(features.TypeIdx(ff)) ;
                    dN_de = face.evalJacobianAt(face.circumcenter) ; % [1 face.nNodes face.nDims(==2)]
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
    
    function elems = slice(this,crossEdg)
    % Slice the element given a signed number of crossed edges
    % signed logical: 0 if the edge is not crossed, or (+/-1) for the sign of the slicing function
        if size(crossEdg,2)~=this.nEdges
            error('Incorrect format for the signed cross edges array') ;
        end
        elems = pkg.mesh.elements.ElementTable ;
        switch this.nDims
            case 1 % slice a bar -> node
                elems = pkg.mesh.elements.ElementTable(...
                                'Types',pkg.mesh.elements.base.Node ...
                                ,'Indices',find(crossEdg(:)).*[0 1] + [1 0] ...
                                ) ;
            case 2 % slice a tri/quad -> bar
                %[elmt,edg,val] = find(crossEdg) ;
                [~,indEdg] = sort(2*abs(crossEdg)+crossEdg,2,'descend') ;
                idx = indEdg(:,1:2) ; % this is a partial solution....
                elems = pkg.mesh.elements.ElementTable(...
                                'Types',pkg.mesh.elements.base.Bar ...
                                ,'Indices',[idx(:,1)*0+1 idx] ...
                                ) ;
            case 3 % slice a 3D elem -> ?
        end
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
        end
    end

    function H = plotReferenceElement(this,ax)
    % Plot the reference element in a given axes
        if nargin<2 ; ax = gca ; end
        h = gobjects(0) ;
        E = this.localCoordinatesIn3D ;
        % Faces
%             h(end+1) = patch(ax,'Vertices',E,'Faces',this.Faces.NodeIdx,'Tag','Faces') ;
%                 h(end).FaceColor = 'w' ;
%                 h(end).EdgeColor = 'k' ;
%                 h(end).LineWidth = 0.5 ;
        h(end+1) = hggroup(ax,'Tag','FaceLabels') ;
            lbl = arrayfun(@(c)['F' num2str(c)],1:this.nFaces,'UniformOutput',false) ;
            P = this.Faces.meanDataAtIndices(E) ;
            txt = text(ax,P(:,1),P(:,2),P(:,3),lbl,'Parent',h(end),'BackgroundColor','w','edgecolor','k') ;
        % Edges
        h(end+1) = patch(ax,'Vertices',E,'Faces',this.Edges.NodeIdx,'Tag','Edges') ;
            h(end).FaceColor = 'none' ;
            h(end).EdgeColor = 'k' ;
            h(end).LineWidth = 2 ;
        h(end+1) = hggroup(ax,'Tag','EdgeLabels') ;
            lbl = arrayfun(@(c)['E' num2str(c)],1:this.nEdges,'UniformOutput',false) ;
            P = this.Edges.meanDataAtIndices(E) ;
            txt = text(ax,P(:,1),P(:,2),P(:,3),lbl,'Parent',h(end),'BackgroundColor','w','edgecolor','k') ;
        % Nodes
        h(end+1) = patch(ax,'Vertices',E,'Faces',(1:this.nNodes)','Tag','Nodes') ;
            h(end).FaceColor = 'none' ;
            h(end).EdgeColor = 'k' ;
            h(end).LineStyle = 'none' ;
            h(end).Marker = '.' ;
            h(end).MarkerSize = 20 ;
        h(end+1) = hggroup(ax,'Tag','NodeLabels') ;
            lbl = arrayfun(@(c)['N' num2str(c)],1:this.nNodes,'UniformOutput',false) ;
            txt = text(ax,E(:,1),E(:,2),E(:,3),lbl,'Parent',h(end),'BackgroundColor','w','edgecolor','k') ;
        % Put in the same group
        H = hggroup(ax) ;
        set(h,'Parent',H) ;
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
                        % Remove the labels
                            delete(findobj(H,'type','text')) ;
                        % Set faces transparent
                            set(findobj(H,'type','patch'),'FaceColor','none')
                        % Highlight the current node only
                            set(findobj(H,'type','patch','-not','marker','none'),'Faces',nn)
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
end

%% ELEMENT IDENTIFICATION FUNCTION
methods
    function bool = isP1Quad(elems)
    % Return true for LagrangeElements('quad',1)
    % Useful to switch node numbering from trigo to P1 quad nodes
        elems = elems(:)' ;
        bool = [elems.nNodes] == 4 ; 
        if any(bool) ; bool = bool & arrayfun(@(e)isa(e,'pkg.mesh.elements.LagrangeElement'),elems) ; end
        if any(bool) ; bool(bool) = bool(bool) & [elems(bool).Order]==1 ; end
        if any(bool) ; bool(bool) = bool(bool) & strcmp({elems(bool).Geometry},'quad') ; end
    end
end

end

