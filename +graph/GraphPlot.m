classdef GraphPlot < handle & matlab.mixin.SetGet & matlab.mixin.Copyable
% Plot a mesh
    
%% PROPERTIES
    properties
        Graph pkg.graph.Graph
        NodeCoordinates
        NodeCoordinatesMode char = 'auto'
        GraphicGroup = gobjects(0)
        ShowLabels = {} ; % label 'edges' and/or 'nodes'
    end
    properties (Access=protected)
        Initialized = false ;
    end
    properties (Dependent)
        % Object
            Parent
            Tag
        % Visibility
            Visible
        % Global Color
            Color
        % Mesh Edges
            EdgeColor
            EdgeWidth
        % Mesh Nodes
            NodeColor
            NodeSize
    end
    properties (Hidden)
        Nodes
        Edges
        Labels
    end
    
%% CONSTRUCTOR/DESTRUCTOR
    methods
        function this = GraphPlot(varargin)
        % Class constructor
            if mod(nargin,2) ; error('Wrong number of arguments') ; end
            this.init ;
            for arg = 1:2:nargin-1
                this.(varargin{arg}) = varargin{arg+1} ;
            end
            this.Initialized = true ;
            this.update ;
        end
        
        function delete(this)
        % Class destructor
            delete(this.GraphicGroup)
        end
    end
    
%% PROPERTIES SET/GET
    methods
    % Changing the mesh ..
        function set.Graph(this,g)
            this.Graph = g ;
            if strcmp(this.NodeCoordinatesMode,'auto')
            % Auto-compute the node coordinates (will update the plot)
                this.NodeCoordinates = this.Graph.autoNodeCoordinates ;
            % re-set the mode to 'auto'
                this.NodeCoordinatesMode = 'auto' ;
            else
                this.update ;
            end
        end
    % Node position
        function set.NodeCoordinates(this,pos)
            this.NodeCoordinates = pos ;
            this.NodeCoordinatesMode = 'manual' ;
            this.update ;
        end
    % Feature labelling
        function set.ShowLabels(this,opt)
            if ischar(opt) ; opt = {opt} ; end
            this.ShowLabels = opt ; 
            update(this) ; 
        end
    end
    
    methods
    % Parent
        function set.Parent(this,parent) ; set(this.GraphicGroup,'Parent',parent) ; end
        function parent = get.Parent(this) ; parent = get(this.GraphicGroup,'Parent') ; end
    % Tag
        function set.Tag(this,tag) ; set(this.GraphicGroup,'Tag',tag) ; end
        function tag = get.Tag(this) ; tag = get(this.GraphicGroup,'Tag') ; end
    % Visibility
        function set.Visible(this,vis) ; set(this.GraphicGroup,'Visible',vis) ; end
        function vis = get.Visible(this) ; vis = get(this.GraphicGroup,'Visible') ; end
    % Global color
        function set.Color(this,clr) ; set(this,'EdgeColor',clr,'NodeColor',clr) ; end
        function clr = get.Color(this) ; clr = get(this,'EdgeColor') ; end
    % Mesh Edges
        function set.EdgeColor(this,clr) ; set(this.Edges,'Color',clr) ; end
        function clr = get.EdgeColor(this) ; clr = get(this.Edges,'Color') ; end
        function set.EdgeWidth(this,w) ; set(this.Edges,'LineWidth',w) ; end
        function w = get.EdgeWidth(this) ; w = get(this.Edges,'LineWidth') ; end
    % Mesh Nodes
        function set.NodeColor(this,clr) ; set(this.Nodes,'MarkerFaceColor',clr,'MarkerEdgeColor',clr) ; end
        function clr = get.NodeColor(this) ; clr = get(this.Nodes,'MarkerEdgeColor') ; end
        function set.NodeSize(this,s) ; set(this.Nodes,'MarkerSize',s) ; end
        function s = get.NodeSize(this) ; s = get(this.Nodes,'MarkerSize') ; end
    end
    
%% GRAPHICS HANDLING
    methods %(Access = protected)
        function init(this)
        % Object initialization function
            this.GraphicGroup = hggroup(gobjects(0)) ;
                this.GraphicGroup.UserData = this ; % Keep a copy of the handle to prevent unwanted deletion
            this.Edges = quiver3([],[],[],[],[],[] ... empty x,y,z,u,v,w data
                            ,0 ... no scaling
                            ,'Color','k' ...
                            ,'LineWidth',1 ...
                            ,'MaxHeadSize',0.075 ...
                            ) ;
            this.Edges.Parent = this.GraphicGroup ;
            this.Nodes = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'Visible','on' ...
                            ,'FaceColor','none' ...
                            ,'EdgeColor','none' ...
                            ,'Marker','o' ...
                            ,'MarkerEdgeColor','k' ...
                            ,'MarkerFaceColor','w' ...
                            ,'MarkerSize',10 ...
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
            if ~this.Initialized ; return ; end
            if isempty(this.Graph) ; return ; end
            % Add dummy nodes for display purposes
                vertices = this.NodeCoordinates ;
                vertices(end+1,:) = NaN ;
                vertices = padarray(vertices,[0 3-size(vertices,2)],0,'post') ;
                edgIdx = [this.Graph.Edges ones(this.Graph.nEdges,1).*this.Graph.nNodes+1] ;
            % Edges
                de = vertices(edgIdx(:,2),:) - vertices(edgIdx(:,1),:) ;
                this.Edges.XData = vertices(edgIdx(:,1),1) ;
                this.Edges.YData = vertices(edgIdx(:,1),2) ;
                this.Edges.ZData = vertices(edgIdx(:,1),3) ;
                this.Edges.UData = de(:,1) ;
                this.Edges.VData = de(:,2) ;
                this.Edges.WData = de(:,3) ;
                this.Edges.ShowArrowHead = this.Graph.Directed ;
            % Nodes
                this.setPatch(this.Nodes,vertices) ;
            % Labels (reinitialized each time...)
                delete(this.Labels)
                this.Labels = gobjects(0) ;
                for ii = 1:numel(this.ShowLabels)
                % Indices, positions, colors
                    switch this.ShowLabels{ii}
                        case 'Nodes'
                            ind = 1:this.Graph.nNodes ;
                            pos = vertices ;
                            color = 'b' ;
                        case 'Edges'
                            ind = 1:this.Graph.nEdges ;
                            pos = 0.5*(vertices(edgIdx(:,1),:) + vertices(edgIdx(:,2),:)) ;
                            color = 'r' ;
                    end
                    if isempty(ind) ; continue ; end
                % 3D position
                    pos = [pos zeros(size(pos,1),3-size(pos,2))] ;
                % Merge superimposed indices
                    [pos,ia,~] = uniquetol(pos(ind,:),'ByRows',true,'OutputAllIndices',true) ;
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
        % Set a patch vertices and faces while reducing strored data
            if nargin<4 % show vertices
                faces = 1:size(vertices,1) ;
            else % reduce the number of vertices by taking only the used ones
                faces(faces>size(vertices,1)) = NaN ;
                [vv,~,ic] = unique(faces(faces>0)) ;
                faces(faces>0) = ic ;
                vertices = vertices(vv,:) ;
            end
            vertices = padarray(vertices,[0 3-size(vertices,2)],0,'post') ;
            set(pa,'vertices',vertices,'faces',faces) ;
        end
    end

    
%% COPY FUNCTION
    methods (Access = protected)
        function obj = copyElement(this)
            obj.GraphicGroup = copy(this.GraphicGroup) ; 
        end
    end


end

