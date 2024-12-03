classdef MeshPlot < handle & matlab.mixin.SetGet & matlab.mixin.Copyable
% Plot a mesh
    
%% PROPERTIES
    properties
        Mesh = pkg.geometry.mesh.Mesh.empty
        NodeCoordinates = [] ;
        GraphicGroup = gobjects(0)
        VisibleNodes = 'none' ;
        VisibleEdges = 'outer' ; 
        VisibleFaces = 'outer' ;
        HighlightBoundaryEdges = true ;
        HighlightEndNodes = false ;
        ShowLabels = {} ;
        ShowFrames = {} ;
        ShowNormals = {} ;
        Selected = struct('Nodes',[],'Elems',[],'Faces',[],'Edges',[]) ;
        Highlighted = struct('Nodes',[],'Elems',[],'Faces',[],'Edges',[]) ;
    end
    properties (Hidden,Access=protected)
        Initialized = false ;
    end
    properties (Dependent)
        % Object
            Parent
            Tag
        % Mesh Deformation (Nodal Displacements)
            Deformation 
        % Visibility
            Visible
        % Global Color
            Color
        % Mesh Faces
            FaceColor
            FaceAlpha
        % Mesh Edges
            EdgeColor
            EdgeWidth
            EdgeAlpha
            EdgeStyle
            EdgeWidthMultiplier
        % Mesh Nodes
            NodeColor
            NodeSize
            NodeStyle
            NodeSizeMultiplier
    end
    properties
    % Color data
        CData
    % Alpha data
        AlphaData
    end
    properties (Hidden)
    % Graphical objects handles
        Nodes
        Edges
        Faces
        EndNodes
        BoundaryEdges
        Frames
        Normals
        Labels
    % Plot vertices
        Vertices
    % Special feature indices
        EdgeIdx % with NaNs at the end)
        FaceTable pkg.geometry.mesh.elements.ElementTable % may be modified to represent simple patches
    % Visible feature indices
        VisibleNodeIdx
        VisibleEdgeIdx
        VisibleFaceIdx
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
            this.Initialized = true ;
            this.update ;
        end
        
        function delete(this)
        % Class destructor
            this.stopSelection ;
            delete(this.GraphicGroup)
        end
    end
    
%% PROPERTIES SET/GET
    methods
    % Changing the mesh ..
        function set.Mesh(this,mesh)
            this.Mesh = mesh ;
            this.NodeCoordinates = mesh.Nodes ;
            this.UpdateOnMeshChange = this.UpdateOnMeshChange ;
            this.update ;
        end
    % Changing the node coordinates
        function set.NodeCoordinates(this,coord)
            if size(coord,1)~=this.Mesh.nNodes ; error('Wrong size') ; end
            coord = coord(:,1:min(end,3)) ;
            coord = padarray(coord,[0 3-size(coord,2)],0,'post') ;
            this.NodeCoordinates = coord ;
            this.update ;
        end
        function coord = get.NodeCoordinates(this)
        % Reset deformation if the mesh has changed !
            if size(this.NodeCoordinates,1)~=this.Mesh.nNodes
                this.NodeCoordinates = this.Mesh.Nodes ;
            end
            coord = this.NodeCoordinates ;
        end
    % Mesh deformation (nodal displacements)
        function u = get.Deformation(this)
            pad = size(this.NodeCoordinates,2)-this.Mesh.nCoord ;
            u = this.NodeCoordinates-padarray(this.Mesh.Nodes,[0 pad],0,'post') ;
        end
        function set.Deformation(this,u)
            pad = size(u,2)-this.Mesh.nCoord ;
            this.NodeCoordinates = padarray(this.Mesh.Nodes,[0 pad],0,'post') + u ;
            if isempty(this.CData) ; this.CData = sqrt(sum(u.^2,2)) ; end
        end
    % Changing the face color
        function cdata = get.CData(this)
            cdata = this.CData ; %this.Faces.FaceVertexCData ;
        end
        function set.CData(this,cdata)
            cdata = full(cdata) ;
            switch size(cdata,1)
                case 0 % no color
                    this.Faces.FaceVertexCData = [] ;
                    this.Edges.FaceVertexCData = [] ;
                    this.CData = [] ; return ;
                case 1 % uniform color, RBG, or string
                    this.Faces.FaceColor = cdata ;
                case this.Mesh.nNodes % nodal values
                    this.Faces.FaceColor = 'interp' ;
                case this.Mesh.nElems % elem values
                    this.Faces.FaceColor = 'flat' ;
                    this.Faces.FaceVertexCData = cdata ;
                    this.CData = [] ; return ;
                case this.Mesh.nEdges % edge values (not working properly...)
                    this.Edges.EdgeColor = 'flat' ;
                    % MATLAB set the edge color to the first node.. 
                    firstnod = this.Mesh.Edges.NodeIdx(:,1) ;
                    clr = NaN(this.Mesh.nNodes+1,size(cdata,2)) ;
                    clr(firstnod,:) = cdata ;
                    this.Edges.FaceVertexCData = clr ;
                    this.CData = [] ; return ;
                otherwise % try on gauss points
                    [ee,~,ie] = this.Mesh.integration ;
                    if size(cdata,1)~=numel(ie) ; error('Wrong shape for CData') ; end
                    this.Faces.FaceColor = 'interp' ;
                    cdata = this.Mesh.interpMat(ee,ie)\cdata ;
            end
            this.CData = cdata ;
            this.update('CData') ;
        end
    % Changing the face transparency
    function adata = get.AlphaData(this)
            adata = this.AlphaData ; %this.Faces.FaceVertexCData ;
        end
        function set.AlphaData(this,adata)
            adata = full(adata) ;
            switch size(adata,1)
                case 0 % no color
                    this.Faces.FaceVertexAlphaData = [] ;
                    this.AlphaData = [] ; return ;
                case 1 % uniform color, RBG, or string
                    this.Faces.FaceAlpha = adata ;
                case this.Mesh.nNodes % nodal values
                    this.Faces.FaceAlpha = 'interp' ;
                case this.Mesh.nElems % elem values
                    this.Faces.FaceAlpha = 'flat' ;
                    this.Faces.FaceVertexAlphaData = adata ;
                    this.AlphaData = [] ; return ;
                otherwise % try on gauss points
                    [ee,~,ie] = this.Mesh.integration ;
                    if size(adata,1)~=numel(ie) ; error('Wrong shape for AlphaData') ; end
                    this.Faces.FaceAlpha = 'interp' ;
                    adata = this.Mesh.interpMat(ee,ie)\adata ;
            end
            this.AlphaData = adata ;
            this.update('AlphaData') ;
        end
    % Feature visibility
        function set.VisibleNodes(this,opt) ; this.VisibleNodes = char(opt) ; update(this) ; end
        function set.VisibleEdges(this,opt) ; this.VisibleEdges = char(opt) ; update(this) ; end
        function set.VisibleFaces(this,opt) ; this.VisibleFaces = char(opt) ; update(this) ; end
        function setVisibleNodeIdx(this)
            switch this.VisibleNodes
                case 'all'
                    this.VisibleNodeIdx = true(this.Mesh.nNodes,1) ;
                case 'outer'
                    this.VisibleNodeIdx = this.Mesh.outerNodes ;
                case 'boundary'
                    this.VisibleNodeIdx = this.Mesh.boundaryNodes ;
                case 'end'
                    this.VisibleNodeIdx = this.Mesh.endNodes ;
                otherwise
                    this.VisibleNodeIdx = false(this.Mesh.nNodes,1) ;
            end
        end
        function setVisibleEdgeIdx(this)
            switch this.VisibleEdges
                case 'all'
                    this.VisibleEdgeIdx = true(this.Mesh.nEdges,1) ;
                case 'outer'
                    this.VisibleEdgeIdx = this.Mesh.outerEdges ;
                case 'boundary'
                    this.VisibleEdgeIdx = this.Mesh.boundaryEdges ;
                otherwise
                    this.VisibleEdgeIdx = false(this.Mesh.nEdges,1) ;
            end
        end
        function setVisibleFaceIdx(this)
        % Empty face list ?
            if isempty(this.Mesh.Faces)
                this.VisibleFaceIdx = false(this.Mesh.nFaces,1) ;
                this.FaceTable = [] ;
                return ;
            end
        % Set visible face indices
            switch this.VisibleFaces
                case 'all'
                    this.VisibleFaceIdx = true(this.Mesh.nFaces,1) ;
                case 'outer'
                    this.VisibleFaceIdx = this.Mesh.outerFaces ;
                otherwise
                    this.VisibleFaceIdx = false(this.Mesh.nFaces,1) ;
            end
        % Compute the element table of visible faces
            this.FaceTable = this.Mesh.Faces.subpart(this.VisibleFaceIdx) ;
            if ~isa(this.Mesh.Elems.Types,'pkg.geometry.mesh.elements.base.BaseElement')
                this.FaceTable = this.FaceTable.simplex ; % to display complicated elements
            end
        end
    % Feature highlighting
        function set.HighlightEndNodes(this,opt) ; this.HighlightEndNodes = logical(opt) ; update(this) ; end
        function set.HighlightBoundaryEdges(this,opt) ; this.HighlightBoundaryEdges = logical(opt) ; update(this) ; end
    % Feature labelling
        function set.ShowLabels(this,opt)
            if ischar(opt) ; opt = {opt} ; end
            this.ShowLabels = opt ; 
            update(this) ; 
        end
    % Local feature frames
        function set.ShowFrames(this,opt)
            if ischar(opt) ; opt = {opt} ; end
            this.ShowFrames = opt ; 
            update(this) ; 
        end
    % Local feature normals
        function set.ShowNormals(this,opt)
            if ischar(opt) ; opt = {opt} ; end
            this.ShowNormals = opt ; 
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
        function set.Color(this,clr) ; set(this,'FaceColor',clr,'EdgeColor',clr,'NodeColor',clr) ; end
        function clr = get.Color(this) ; clr = get(this,'EdgeColor') ; end
    % Mesh Faces
        function set.FaceColor(this,clr) ; set(this.Faces,'FaceColor',clr) ; end
        function clr = get.FaceColor(this) ; clr = get(this.Faces,'FaceColor') ; end
        function set.FaceAlpha(this,a) ; set(this.Faces,'FaceAlpha',a) ; end
        function a = get.FaceAlpha(this) ; a = get(this.Faces,'FaceAlpha') ; end
    % Mesh Edges
        function set.EdgeStyle(this,sty) ; set([this.Edges this.BoundaryEdges],'LineStyle',sty) ; end
        function sty = get.EdgeStyle(this) ; sty = get(this.Edges,'LineStyle') ; end
        function set.EdgeColor(this,clr) ; set([this.Edges this.BoundaryEdges],'EdgeColor',clr) ; end
        function clr = get.EdgeColor(this) ; clr = get(this.Edges,'EdgeColor') ; end
        function set.EdgeAlpha(this,a) ; set([this.Edges this.BoundaryEdges],'EdgeAlpha',a) ; end
        function a = get.EdgeAlpha(this) ; a = get(this.Edges,'EdgeAlpha') ; end
        function set.EdgeWidth(this,w) ; set([this.Edges this.BoundaryEdges],{'LineWidth'},{w ; w*this.EdgeWidthMultiplier}) ; end
        function w = get.EdgeWidth(this) ; w = get(this.Edges,'LineWidth') ; end
        function set.EdgeWidthMultiplier(this,m) ; w = this.EdgeWidth ; set([this.Edges this.BoundaryEdges],{'LineWidth'},{w ; w*m}) ; end
        function m = get.EdgeWidthMultiplier(this) ; m = get(this.BoundaryEdges,'LineWidth')/get(this.Edges,'LineWidth') ; end
    % Mesh Nodes
        function set.NodeColor(this,clr) ; set([this.Nodes this.EndNodes],'MarkerFaceColor',clr,'MarkerEdgeColor',clr) ; end
        function clr = get.NodeColor(this) ; clr = get(this.Nodes,'MarkerEdgeColor') ; end
        function set.NodeSize(this,s) ; set([this.Nodes this.EndNodes],{'MarkerSize'},{s ; s*this.NodeSizeMultiplier}) ; end
        function s = get.NodeSize(this) ; s = get(this.Nodes,'MarkerSize') ; end
        function set.NodeSizeMultiplier(this,m) ; s = this.NodeSize ; set([this.Nodes this.EndNodes],{'MarkerSize'},{s ; s*m}) ; end
        function m = get.NodeSizeMultiplier(this) ; m = get(this.EndNodes,'MarkerSize')/get(this.Nodes,'MarkerSize') ; end
    end
    
%% UPDATE ON MESH CHANGES
    properties
        UpdateOnMeshChange = false
    end
    properties (Hidden,Transient)
        MeshUpdateListeners = event.listener.empty
    end
    methods
        function set.UpdateOnMeshChange(this,bool)
            delete(this.MeshUpdateListeners)
            if bool % re-start update on mesh changes
                this.MeshUpdateListeners(end+1) = addlistener(this.Mesh,'Nodes','PostSet',@this.meshUpdateCallback) ;
                this.MeshUpdateListeners(end+1) = addlistener(this.Mesh,'Elems','PostSet',@this.meshUpdateCallback) ;
            end
            this.UpdateOnMeshChange = bool ;
        end
        function meshUpdateCallback(this,varargin)
            this.update ;
            %drawnow ;
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
            this.Frames = gobjects(0) ;
            this.Normals = gobjects(0) ;
            this.Labels = text(gobjects(0)...
                            ,NaN,NaN,'' ...
                            ,'Parent',this.GraphicGroup ...
                            ,'Visible','off' ...
                            ) ;
        end
        
        function update(this,what)
        % Object updating function
            if ~this.Initialized ; return ; end
            if isempty(this.Mesh) ; return ; end
            if nargin<2 ; what = 'all' ; end
            if isempty(what) ; return ; end
            toUpdate = @(arg)any(ismember(arg,what)) ;
            % Add dummy nodes for display purposes
                if toUpdate({'all','Nodes','Vertices'})
                    this.Vertices = this.NodeCoordinates ;
                    this.Vertices(end+1,:) = NaN ;
                end
            % Faces
                if toUpdate({'all','Faces'}) 
                    this.setVisibleFaceIdx ; 
                end
                if toUpdate({'all','Faces','Vertices','CData','AlphaData'}) 
                    this.setPatch(this.Faces,this.Vertices,this.FaceTable.indicesWithNaNs,this.CData,this.AlphaData) ;
                end
            % Edges
                if toUpdate({'all','Edges'}) 
                    this.EdgeIdx = [this.Mesh.Edges.indicesWithNaNs ones(this.Mesh.nEdges,1).*this.Mesh.nNodes+1] ;
                    this.setVisibleEdgeIdx ;
                end
                if toUpdate({'all','Edges','Vertices'}) 
                    this.setPatch(this.Edges,this.Vertices,this.EdgeIdx(this.VisibleEdgeIdx,:)) ;
                end
            % Highlighted Edges
                if toUpdate({'all','Edges'}) 
                    if this.HighlightBoundaryEdges
                        highlightedEdges = this.Mesh.boundaryEdges ;
                    else
                        highlightedEdges = [] ;
                    end
                end
                if toUpdate({'all','Edges','Vertices'}) 
                    this.setPatch(this.BoundaryEdges,this.Vertices,this.EdgeIdx(highlightedEdges,:)) ;
                end
            % Nodes
                if toUpdate({'all','Nodes'}) 
                    this.setVisibleNodeIdx ;
                end
                if toUpdate({'all','Nodes','Vertices'}) 
                    this.setPatch(this.Nodes,this.Vertices(this.VisibleNodeIdx,:)) ;
                end
            % Highlighted Nodes
                if toUpdate({'all','Nodes'})
                    if this.HighlightEndNodes
                        highlightedNodes = reshape(find(this.Mesh.endNodes),1,[]) ;
                    else
                        highlightedNodes = [] ;
                    end
                end
                if toUpdate({'all','Nodes','Vertices'})
                    this.setPatch(this.EndNodes,this.Vertices(highlightedNodes,:)) ;
                end
            % Frames
                if toUpdate({'all','Frames','Vertices'})
                    delete(this.Frames)
                    this.Frames = gobjects(0) ;
                    % Plot the frames
                        for ii = 1:numel(this.ShowFrames)
                            switch this.ShowFrames{ii}
                                case 'Nodes'
                                    origin = padarray(this.Mesh.Nodes,[0 3*this.Mesh.nCoord],0,'post') ;
                                    frames = this.Mesh.nodeNormals ;
                                otherwise
                                    origin = padarray(this.Mesh.centroid(this.ShowFrames{ii}),[0 3*this.Mesh.nCoord],0,'post') ;
                                    frames = permute(this.Mesh.getFrames(this.ShowFrames{ii}),[3 2 1]) ;
                            end
                            nVec = size(frames,3) ;
                            colors = linspecer(nVec) ;
                            for vv = 1:nVec
                                this.Frames(end+1) = quiver3( ...
                                                        origin(:,1),origin(:,2),origin(:,3) ...
                                                        ,frames(:,1,vv),frames(:,2,vv),frames(:,3,vv) ...
                                                        ,'Color',colors(vv,:) ...
                                                        ,'Parent',this.GraphicGroup ...
                                                        ) ;
                            end
                        end
                    % Set common properties
                        set(this.Frames ...
                            ,'LineWidth',1 ...
                            ,'MaxHeadSize',0.05 ...
                            ,'AutoScaleFactor',0.15 ...
                            ,'Parent',this.GraphicGroup ...
                            ) ;
                end
            % Normals
                if toUpdate({'all','Normals','Vertices'})
                    delete(this.Normals)
                    this.Normals = gobjects(0) ;
                    % Plot the frames
                        for ii = 1:numel(this.ShowNormals)
                            switch this.ShowNormals{ii}
                                case 'Nodes'
                                    origin = padarray(this.Mesh.Nodes,[0 3-this.Mesh.nCoord],0,'post') ;
                                    normals = this.Mesh.nodeNormals ;
                                case 'BoundaryNodes'
                                    [~,~,normals,no] = this.Mesh.boundaryNormals ;
                                    origin = padarray(this.Mesh.Nodes(no,:),[0 3-this.Mesh.nCoord],0,'post') ;
                                case 'BoundaryEdges'
                                    [normals,ed] = this.Mesh.boundaryNormals ;
                                    origin = padarray(this.Mesh.centroid(this.Mesh.Edges.subpart(ed)),[0 3-this.Mesh.nCoord],0,'post') ;
                                otherwise
                                    origin = padarray(this.Mesh.centroid(this.ShowNormals{ii}),[0 3-this.Mesh.nCoord],0,'post') ;
                                    normals = this.Mesh.getNormals(this.ShowNormals{ii}) ;
                            end
                            this.Normals(end+1) = quiver3( ...
                                                    origin(:,1),origin(:,2),origin(:,3) ...
                                                    ,normals(:,1),normals(:,2),normals(:,3) ...
                                                    ) ;
                        end
                    % Set common properties
                        set(this.Normals ...
                            ,'Color','k' ...
                            ,'LineWidth',1 ...
                            ,'MaxHeadSize',0.05 ...
                            ,'AutoScaleFactor',0.15 ...
                            ,'Parent',this.GraphicGroup ...
                            ) ;
                end
            % Labels (reinitialized each time...)
                if toUpdate({'all','Labels','Vertices'})
                    delete(this.Labels)
                    this.Labels = gobjects(0) ;
                    for ii = 1:numel(this.ShowLabels)
                    % Indices, positions, colors
                        switch this.ShowLabels{ii}
                            case 'Nodes'
                                ind = find(this.VisibleNodeIdx) ;
                                pos = this.Mesh.Nodes ;
                                color = 'b' ;
                            case 'Edges'
                                ind = find(this.VisibleEdgeIdx) ;
                                pos = this.Mesh.centroid(this.Mesh.Edges) ;
                                color = 'r' ;
                            case 'Faces'
                                ind = find(this.VisibleFaceIdx) ;
                                pos = this.Mesh.centroid(this.Mesh.Faces) ;
                                color = 'm' ;
                            case 'Elems'
                                ind = 1:this.Mesh.nElems ;
                                pos = this.Mesh.centroid(this.Mesh.Elems) ;
                                color = 'k' ;
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
    end
    
    
%% UNIQUE INDICES (DATA REDUCTION)
    methods
        function setPatch(~,pa,vertices,faces,cdata,adata)
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
        % Color data
            if nargin<5 || isempty(cdata) ; return ; end
            set(pa,'facevertexcdata',cdata(vv,:)) ;
        % Alpha data
            if nargin<6 || isempty(adata) ; return ; end
            set(pa,'facevertexalphadata',adata(vv,:)) ;
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
            vertices = padarray(this.NodeCoordinates,[1 0],NaN,'post') ; 
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
            vertices = padarray(this.NodeCoordinates,[1 0],NaN,'post') ; 
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
                                dist = sqrt(sum((this.Mesh.Nodes-pos).^2,2)) ;
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