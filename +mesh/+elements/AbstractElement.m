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
        
%         function [D,ORD] = evalAllDerivativesAt(this,E,T,delta)
%         % Evaluate all derivatives corresponding to the total order T
%         % E = [nE nDims]
%         % delta = optionnal, for finite diff derivatives (see above)
%         % D = [nE nNodes nDeriv]
%         % So that d^{T}N(E)_(de_i1.de_i2...de_iT) = D(:,:,all(ORD==[i1,i2,...,iT],2))*Ne
%             if nargin<3 ; T = 1 ; end % Return the gradient by default
%             if nargin<4 ; delta = this.defaultDelta ; end
%         % Build the derivation orders
%             ORD = perms(repmat(0:T,[1 this.nDims])) ;
%             ORD = unique(ORD(:,1:this.nDims),'rows') ;
%             ORD = ORD(sum(ORD,2)==T,:) ;
%             [~,sortInd] = sort(sum((ORD+1).^(this.nDims:-1:1),2),'descend') ; 
%             ORD = ORD(sortInd,:) ;
%             nDeriv = size(ORD,1) ;
%         % Initialize the matrix
%             D = NaN([size(E,1) this.nNodes nDeriv]) ;
%         % Compute
%             for oo = 1:nDeriv
%                 D(:,:,oo) = this.evalDerivativeAt(E,ORD(oo,:),delta) ;
%             end
%         end
        
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
            subdiv = 100*ones(1,this.nDims) ; % Number of Subdivisions
            switch this.nDims % Displayed derivation orders
                case 1
                    ORD = [0;1;2] ;
                case 2
                    ORD = [0 0 ; 1 0 ; 0 1 ; 2 0 ; 0 2 ; 1 1] ;
                case 3
            end
            % Local coordinates
                bbox = this.localCoordinatesDomain ;
                E = arrayfun(@(dim)linspace(bbox(1,dim),bbox(2,dim),subdiv(dim)),1:this.nDims,'UniformOutput',false) ;
                if numel(E)>1 ; [E{:}] = ndgrid(E{:}) ; end
                EE = cat(this.nDims+1,E{:}) ;
            % Figure
                clf(fig) ;
                H = gobjects(0) ;
                ax = gobjects(0) ;
                lin = gobjects(0) ;
                srf = gobjects(0) ;
                ttl = gobjects(0) ;
             % For each derivation order
                nOrd = size(ORD,1) ;
                for oo = 1:nOrd
                    % Evaluate the shape function derivatives
                        e = reshape(EE,[],this.nDims) ;
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
                            % Highlight the current point
                                switch this.nDims
                                    case 1
                                        lin(end+1) = plot(ax(end),EE(:,:,1),NN(:,nn)) ;
                                    case 2
                                        srf(end+1) = surf(ax(end),EE(:,:,1),EE(:,:,2),NN(:,:,nn),NN(:,:,nn)) ;
                                    case 3
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
    
end

