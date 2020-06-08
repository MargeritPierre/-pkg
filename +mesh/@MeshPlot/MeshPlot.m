classdef MeshPlot < handle & matlab.mixin.SetGet & matlab.mixin.Copyable
% Plot a mesh
    
%% PROPERTIES
    properties
        Mesh = pkg.mesh.Mesh.empty
        CData = [] 
        GraphicGroup = gobjects(0)
        VisibleNodes = 'none' ;
        VisibleEdges = 'outer' ; 
        VisibleFaces = 'outer' ;
        HighlightBoundaryEdges = true ;
        HighlightEndNodes = false ;
        ShowLabels = {} ;
        Selected = struct('Nodes',[],'Elems',[],'Faces',[],'Edges',[]) ;
        Highlighted = struct('Nodes',[],'Elems',[],'Faces',[],'Edges',[]) ;
    end
    properties (Dependent)
        % Object
            Parent
            Tag
        % Mesh Faces
            FaceColor
            FaceAlpha
        % Mesh Edges
            EdgeColor
            EdgeWidth
            EdgeAlpha
            EdgeStyle
        % Mesh Nodes
            NodeColor
            NodeSize
            NodeStyle
    end
    properties (Hidden)
        Nodes
        Edges
        Faces
        EndNodes
        BoundaryEdges
        Labels
    end
    
%% CONSTRUCTOR/DESTRUCTOR
    methods
        function this = MeshPlot(varargin)
        % Class constructor
            if mod(nargin,2) ; error('Wrong number of arguments') ; end
            this.init ;
            for arg = 1:2:nargin-1
                this.(varargin{arg}) = varargin{arg+1} ;
            end
        end
        
        function delete(this)
        % Class destructor
            this.stopSelection ;
            delete(this.GraphicGroup)
        end
    end
    
%% PROPERTIES SET/GET
    methods
        function set.Mesh(this,mesh)
        % On the mesh set
            this.Mesh = mesh ;
            this.update ;
        end
        function set.CData(this,cdata)
        end
        function set.VisibleNodes(this,opt) ; this.VisibleNodes = char(opt) ; update(this) ; end
        function set.VisibleEdges(this,opt) ; this.VisibleEdges = char(opt) ; update(this) ; end
        function set.VisibleFaces(this,opt) ; this.VisibleFaces = char(opt) ; update(this) ; end
        function set.HighlightEndNodes(this,opt) ; this.HighlightEndNodes = logical(opt) ; update(this) ; end
        function set.HighlightBoundaryEdges(this,opt) ; this.HighlightBoundaryEdges = logical(opt) ; update(this) ; end
        function set.ShowLabels(this,opt)
            if ischar(opt) opt = {opt} ; end
            this.ShowLabels = opt ; 
            update(this) ; 
        end
    end
    
    methods
        function set.Parent(this,parent)
            set(this.GraphicGroup,'Parent',parent) ;
        end
        function parent = get.Parent(this)
            parent = get(this.GraphicGroup,'Parent') ;
        end
    end
    
    
%% GRAPHICS HANDLING
    methods %(Access = protected)
        function init(this)
        % Object initialization function
            this.GraphicGroup = hggroup(gobjects(0)) ;
                this.GraphicGroup.UserData = this ; % Keep a copy of the handle to prevent unwanted deletion
            this.Faces = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'Visible','on' ...
                            ,'FaceColor','w' ...
                            ,'EdgeColor','none' ...
                            ,'Parent',this.GraphicGroup ...
                            ) ;
            this.Edges = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'FaceColor','none' ...
                            ,'EdgeColor','k' ...
                            ,'Parent',this.GraphicGroup ...
                            ) ;
            this.BoundaryEdges = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'FaceColor','none' ...
                            ,'EdgeColor','k' ...
                            ,'LineWidth',1.5 ...
                            ,'Parent',this.GraphicGroup ...
                            ) ;
            this.Nodes = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'Visible','on' ...
                            ,'FaceColor','none' ...
                            ,'EdgeColor','none' ...
                            ,'Marker','.' ...
                            ,'MarkerEdgeColor','k' ...
                            ,'MarkerFaceColor','k' ...
                            ,'MarkerSize',6 ...
                            ,'Parent',this.GraphicGroup ...
                            ) ;
            this.EndNodes = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'Visible','on' ...
                            ,'FaceColor','none' ...
                            ,'EdgeColor','none' ...
                            ,'Marker','.' ...
                            ,'MarkerEdgeColor','k' ...
                            ,'MarkerFaceColor','k' ...
                            ,'MarkerSize',15 ...
                            ,'Parent',this.GraphicGroup ...
                            ) ;
            this.Labels = text(gobjects(0)...
                            ,NaN,NaN,'' ...
                            ,'Parent',this.GraphicGroup ...
                            ,'Visible','off' ...
                            ) ;
        end
        
        function update(this)
        % Object updating function
            % Add dummy nodes for display purposes
                vertices = [this.Mesh.X.Values ; NaN(1,this.Mesh.nCoord)] ;
                edgIdx = [this.Mesh.Edges.indicesWithNaNs ones(this.Mesh.nEdges,1)*this.Mesh.nNodes+1] ;
            % Set vertices to patches
                %ppp = findobj(this.GraphicGroup,'type','patch') ;
                %set(ppp,'Vertices',vertices) ;
            % Faces
                switch this.VisibleFaces
                    case 'all'
                        visibleFaces = true(this.Mesh.nFaces,1) ;
                    case 'outer'
                        visibleFaces = this.Mesh.outerFaces ;
                    otherwise
                        visibleFaces = false(this.Mesh.nFaces,1) ;
                end
                this.setPatch(this.Faces,vertices,this.Mesh.Faces.NodeIdx(visibleFaces,:)) ;
            % Edges
                switch this.VisibleEdges
                    case 'all'
                        visibleEdges = true(this.Mesh.nEdges,1) ;
                    case 'outer'
                        visibleEdges = this.Mesh.outerEdges ;
                    case 'boundary'
                        visibleEdges = this.Mesh.boundaryEdges ;
                    otherwise
                        visibleEdges = false(this.Mesh.nEdges,1) ;
                end
                this.setPatch(this.Edges,vertices,edgIdx(visibleEdges,:)) ;
            % Highlighted Edges
                if this.HighlightBoundaryEdges
                    highlightedEdges = this.Mesh.boundaryEdges ;
                else
                    highlightedEdges = [] ;
                end
                this.setPatch(this.BoundaryEdges,vertices,edgIdx(highlightedEdges,:)) ;
            % Nodes
                switch this.VisibleNodes
                    case 'all'
                        visibleNodes = true(this.Mesh.nNodes,1) ;
                    case 'outer'
                        visibleNodes = this.Mesh.outerNodes ;
                    case 'boundary'
                        visibleNodes = this.Mesh.boundaryNodes ;
                    case 'end'
                        visibleNodes = this.Mesh.endNodes ;
                    otherwise
                        visibleNodes = false(this.Mesh.nNodes,1) ;
                end
                this.setPatch(this.Nodes,vertices(visibleNodes,:)) ;
            % Highlighted Nodes
                if this.HighlightEndNodes
                    highlightedNodes = reshape(find(this.Mesh.endNodes),1,[]) ;
                else
                    highlightedNodes = [] ;
                end
                this.setPatch(this.EndNodes,vertices(highlightedNodes,:)) ;
            % Labels (reinitialized each time...)
                delete(this.Labels)
                this.Labels = gobjects(0) ;
                for ii = 1:numel(this.ShowLabels)
                % Indices, positions, colors
                    switch this.ShowLabels{ii}
                        case 'Nodes'
                            ind = find(visibleNodes) ;
                            pos = this.Mesh.X.Values ;
                            color = 'b' ;
                        case 'Edges'
                            ind = find(visibleEdges) ;
                            pos = this.Mesh.circumCenters(this.Mesh.Edges) ;
                            color = 'r' ;
                        case 'Faces'
                            ind = find(visibleFaces) ;
                            pos = this.Mesh.circumCenters(this.Mesh.Faces) ;
                            color = 'm' ;
                        case 'Elems'
                            ind = 1:this.Mesh.nElems ;
                            pos = this.Mesh.circumCenters(this.Mesh.Elems) ;
                            color = 'k' ;
                    end
                    if isempty(ind) ; continue ; end
                % 3D position
                    pos = [pos zeros(size(pos,1),3-size(pos,2))] ;
                % Merge superimposed indices
                    [pos,ia,ic] = uniquetol(pos(ind,:),'ByRows',true,'OutputAllIndices',true) ;
                    ind = cellfun(@(ii)ind(ii),ia,'UniformOutput',false) ;
                % Label text
                    txt = cellfun(@num2str,ind,'UniformOutput',false) ;
                % Create label
                    lbl = text(pos(:,1),pos(:,2),pos(:,3),txt) ;
                % Set specific properties
                    set(lbl,'Color',color) ;
                % Add to the list
                    this.Labels = [this.Labels ; lbl(:)] ;
                end
            % Set common properties
                set(this.Labels,'Parent',this.GraphicGroup) ;
                set(this.Labels,'FontName','consolas','interpreter','tex') ;
                set(this.Labels,'FontSize',12,'FontWeight','bold') ;
        end
    end
    
    
%% UNIQUE INDICES (DATA REDUCTION)
    methods
        function setPatch(~,pa,vertices,faces)
        % Set a patch vertices and faces while reducing stroed data
            if nargin<4 % show vertices
                faces = 1:size(vertices,1) ;
            else % reduce the number of vertices by taking only the used ones
                [vv,~,ic] = unique(faces(faces>0)) ;
                faces(faces>0) = ic ;
                vertices = vertices(vv,:) ;
            end
            set(pa,'vertices',vertices,'faces',faces) ;
        end
    end

    
%% COPY FUNCTION
    methods (Access = protected)
        function obj = copyElement(this)
            obj.GraphicGroup = copy(this.GraphicGroup) ; 
        end
    end
    
    
%% SELECTION ABILITIES
    properties (Hidden)
        SelectionMode = '' ;
        SelectionGroup = gobjects(0)
        HighlightedNodes = gobjects(0)
        HighlightedFeatures = gobjects(0)
        SelectedNodes = gobjects(0)
        SelectedFeatures = gobjects(0)
        SelectionListeners
        StoredCallback
    end
    methods
        function set.Selected(this,selected)
        % Change the selection status
            if isempty(this.SelectedNodes) ...
                    || ~isvalid(this.SelectedNodes) ...
                this.initSelectionGraphics ;
            end
        % Set the property
            this.Selected = selected ;
        % Nodes
            vertices = [this.Mesh.X.Values ; NaN(1,this.Mesh.nCoord)] ;
            this.setPatch(this.SelectedNodes,vertices(this.Selected.Nodes,:)) ;
        % Other features (Elems, Faces, Edges)
            IDX = [] ;
            for feat = {'Elems' 'Faces' 'Edges'}
                if isempty(this.Selected.(feat{:})) ; continue ; end
                % Indices
                    idx = this.Mesh.(feat{:}).indicesWithNaNs ;
                    idx = idx(this.Selected.(feat{:}),:) ; 
                % Add dummy node to edges
                    if strcmp(feat{:},'Edges') ; idx = [ones(size(idx,1),1)*this.Mesh.nNodes+1 idx] ; end
                % Add to the global list of indices
                    IDX = [IDX NaN(size(IDX,1),size(idx,2)-size(IDX,2))] ;
                    IDX = [IDX ; idx NaN(size(idx,1),size(IDX,2)-size(idx,2))] ;
            end
            this.setPatch(this.SelectedFeatures,vertices,IDX) ;
        end
        
        function set.Highlighted(this,highlighted)
        % Change the selection status
            if isempty(this.HighlightedNodes) ...
                    || ~isvalid(this.HighlightedNodes) ...
                this.initSelectionGraphics ;
            end
        % Set the property
            this.Highlighted = highlighted ;
        % Nodes
            vertices = [this.Mesh.X.Values ; NaN(1,this.Mesh.nCoord)] ;
            this.setPatch(this.HighlightedNodes,vertices(this.Highlighted.Nodes,:)) ;
        % Other features (Elems, Faces, Edges)
            IDX = [] ;
            for feat = {'Elems' 'Faces' 'Edges'}
                if isempty(this.Highlighted.(feat{:})) ; continue ; end
                % Indices
                    idx = this.Mesh.(feat{:}).indicesWithNaNs ;
                    idx = idx(this.Highlighted.(feat{:}),:) ; 
                % Add dummy node to edges
                    if strcmp(feat{:},'Edges') ; idx = [ones(size(idx,1),1)*this.Mesh.nNodes+1 idx] ; end
                % Add to the global list of indices
                    IDX = [IDX NaN(size(IDX,1),size(idx,2)-size(IDX,2))] ;
                    IDX = [IDX ; idx NaN(size(idx,1),size(IDX,2)-size(idx,2))] ;
            end
            this.setPatch(this.HighlightedFeatures,vertices,IDX) ;
        end
        
        function initSelectionGraphics(this)
        % Init the graphical objects used for interactive selection
            this.stopSelection ;
            this.SelectionGroup = hggroup(this.GraphicGroup) ;
            this.SelectedNodes = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'Visible','on' ...
                            ,'FaceColor','none' ...
                            ,'EdgeColor','none' ...
                            ,'Marker','.' ...
                            ,'MarkerEdgeColor','r' ...
                            ,'MarkerFaceColor','r' ...
                            ,'MarkerSize',12 ...
                            ,'Parent',this.SelectionGroup ...
                            ) ;
            this.SelectedFeatures = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'FaceColor','none' ...
                            ,'EdgeColor','r' ...
                            ,'LineWidth',1.5 ...
                            ,'Parent',this.SelectionGroup ...
                            ) ;
            this.HighlightedNodes = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'Visible','on' ...
                            ,'FaceColor','none' ...
                            ,'EdgeColor','none' ...
                            ,'Marker','.' ...
                            ,'MarkerEdgeColor','b' ...
                            ,'MarkerFaceColor','b' ...
                            ,'MarkerSize',12 ...
                            ,'Parent',this.SelectionGroup ...
                            ) ;
            this.HighlightedFeatures = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'FaceColor','none' ...
                            ,'EdgeColor','b' ...
                            ,'LineWidth',1.5 ...
                            ,'Parent',this.SelectionGroup ...
                            ) ;
        end
        
        function select(this,features)
        % Start the selection of features (Nodes, Elems, Faces or Edges)
            if nargin<2 ; features = 'Elems' ; end
            if isempty(features) ; this.stopSelection ; return ; end
        % Change the selection mode
            this.SelectionMode = features ;
        % Delete listeners
            delete(this.SelectionListeners) ;
            this.SelectionListeners = event.listener.empty ;
        % Force intialization of selection graphics
            this.Selected = this.Selected ;
        % Reset listeners
            fig = ancestor(this.GraphicGroup,'figure') ;
            %this.SelectionListeners(end+1) = listener(fig,'WindowMouseMotion',@this.selectionCallback) ;
            this.StoredCallback = fig.WindowButtonMotionFcn ; % we need to do this for the CurrentPoint axe property to be updated....
            fig.WindowButtonMotionFcn = @this.selectionCallback ;
            ax = ancestor(this.GraphicGroup,'axes') ;
            ax.Interactions = ax.Interactions ;
            this.SelectionListeners(end+1) = listener(fig,'WindowMousePress',@this.selectionCallback) ;
        end
        
        function stopSelection(this)
        % Stop the selection of features
            delete(this.SelectionListeners) ;
            delete(this.SelectionGroup) ;
            set(ancestor(this.GraphicGroup,'figure'),'WindowButtonMotionFcn',this.StoredCallback) ;
        end
        
        function clearSelection(this,features)
        % Clear the selected features (Nodes, Elems, Faces or Edges)
            if nargin<2 ; features = fieldnames(this.Selected) ; end
            if ischar(features) ; features = {features} ; end
            for ff = 1:numel(features)
                this.Selected.(features{ff}) = [] ;
                this.Highlighted.(features{ff}) = [] ;
            end
        end
        
        function selectionCallback(this,fig,evt)
        % Responding to a selection listener...
            switch evt.EventName
                case 'WindowMouseMotion'
                    % Get the hitted component
                        ax = ancestor(this.GraphicGroup,'axes') ;
                        pos = ax.CurrentPoint(1,1:2) ;
                        switch this.SelectionMode
                            case 'Nodes'
                                dist = sqrt(sum((this.Mesh.X.Values-pos).^2,2)) ;
                                ind = find(dist<norm(range([ax.XLim ; ax.YLim],2))*5e-3) ;
                            otherwise
                                ind = find(this.Mesh.isInside(pos,this.Mesh.(this.SelectionMode))) ;
                        end
                    % Highlight it
                        this.Highlighted.(this.SelectionMode) = ind ;
                case 'WindowMousePress'
                    feat = this.SelectionMode ;
                    switch fig.SelectionType
                        case 'normal' % Left-click
                            this.Selected.(feat) = this.Highlighted.(feat) ;
                        case 'alt' % Ctrl + Left-click
                            this.Selected.(feat) = setdiff(this.Selected.(feat),this.Highlighted.(feat)) ;
                        case 'extend' % Shift + Left-click
                            this.Selected.(feat) = unique([this.Selected.(feat) this.Highlighted.(feat)]) ;
                    end
            end
        end
    end




end

