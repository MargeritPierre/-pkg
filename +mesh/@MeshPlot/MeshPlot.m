classdef MeshPlot < handle & matlab.mixin.SetGet & matlab.mixin.Copyable
% Plot a mesh
    
%% PROPERTIES
    properties
        Mesh = pkg.mesh.Mesh.empty
        CData = [] 
        GraphicGroup = gobjects(0)
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
        AllEdges
        BoundaryEdges
        OuterFaces
        Labels
    end
    
%% CONSTRUCTOR/DESTRUCTOR
    methods
        function this = MeshPlot(varargin)
        % Class constructor
            if mod(nargin,2) ; error('wrong number of arguments') ; end
            this.init ;
            for arg = 1:2:nargin-1
                this.(varargin{arg}) = varargin{arg+1} ;
            end
        end
        
        function delete(this)
        % Class destructor
%             delete(this.Patches)
%             delete(this.Lines)
%             delete(this.Labels)
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
    end
    
    methods
        function set.Parent(this,parent)
            set(this.GraphicGroup,'Parent',parent) ;
        end
        function parent = get.Parent(this)
            parent = get(this.GraphicGroup,'Parent') ;
        end
    end
    
    methods (Access = protected)
        function setAllProperties(this,varargin)
            %objects = [this.Patch this.
        end
    end
    
    
%% GRAPHICS HANDLING
    methods (Access = protected)
        function init(this)
        % Object initialization function
            this.GraphicGroup = hggroup(gobjects(0)) ;
            this.Nodes = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'Visible','off' ...
                            ,'FaceColor','none' ...
                            ,'EdgeColor','none' ...
                            ,'Marker','.' ...
                            ,'MarkerEdgeColor','k' ...
                            ,'MarkerFaceColor','k' ...
                            ,'MarkerSize',6 ...
                            ,'Parent',this.GraphicGroup ...
                            ) ;
            this.AllEdges = patch(gobjects(0) ...
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
            this.OuterFaces = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'FaceColor','w' ...
                            ,'EdgeColor','none' ...
                            ,'Parent',this.GraphicGroup ...
                            ) ;
            this.Labels = text(gobjects(0)...
                            ,NaN,NaN,'' ...
                            ,'Parent',this.GraphicGroup ...
                            ) ;
        end
        
        function update(this)
        % Object updating function
            % Add dummy nodes for display purposes
                vertices = [this.Mesh.X.Values ; NaN(1,this.Mesh.nCoord)] ;
                edgIdx = [this.Mesh.Edges.indicesWithNaNs ones(this.Mesh.nEdges,1)*this.Mesh.nNodes+1] ;
            % Node plot
                this.Nodes.Vertices = vertices ;
                this.Nodes.Faces = 1:size(vertices,1) ;
            % All Edges
                this.AllEdges.Vertices = vertices ;
                this.AllEdges.Faces = edgIdx ;
            % Boundary Edges
                this.BoundaryEdges.Vertices = vertices ;
                this.BoundaryEdges.Faces = edgIdx(this.Mesh.boundaryEdges,:) ;
        end
    end
    
    
%% COPY FUNCTION
    methods (Access = protected)
        function obj = copyElement(this)
            obj.GraphicGroup = copy(this.GraphicGroup) ; 
        end
    end




end

