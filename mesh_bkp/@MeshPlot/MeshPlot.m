classdef MeshPlot < handle & matlab.mixin.SetGet & matlab.mixin.Copyable
% Plot a mesh
    
%% PROPERTIES
    properties
        Mesh = pkg.mesh.Mesh.empty
        CData = [] 
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
        Patches
        Lines
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
            delete(this.Patches)
            delete(this.Lines)
            delete(this.Labels)
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
            set([this.Patches this.Lines this.Labels],'Parent',parent) ;
        end
        function parent = get.Parent(this)
            parent = get(this.Patches,'Parent') ;
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
            this.Patches = patch(gobjects(0) ...
                            ,'Vertices',[] ...
                            ,'Faces',[] ...
                            ,'FaceColor','none' ...
                            ,'EdgeColor','k' ...
                            ) ;
            this.Lines = line(gobjects(0) ...
                            ,NaN ...
                            ,NaN ...
                            ,'Linewidth',1.5 ...
                            ,'Color','k' ...
                            ) ;
            this.Labels = text(gobjects(0),NaN,NaN,'') ;
        end
        
        function update(this)
        % Object updating function
            % Patch
                this.Patches.Vertices = this.Mesh.Nodes ;
                this.Patches.Faces = this.Mesh.Faces ;
            % Boundary curves
                crv = this.Mesh.boundaryCurves ;
                nans = isnan(crv) ;
                crv(nans) = 1 ;
                crv = this.Mesh.Nodes(crv,:) ;
                crv(nans,:) = NaN ;
                this.Lines.XData = crv(:,1) ;
                this.Lines.YData = crv(:,2) ;
                if size(crv,2)>2 ; this.Lines.ZData = crv(:,3) ; end
        end
    end
    
    
%% COPY FUNCTION
    methods (Access = protected)
        function obj = copyElement(this)
            obj.Patches = copy(this.Patches) ; 
            obj.Lines = copy(this.Lines) ; 
            obj.Labels = copy(this.Labels) ; 
        end
    end




end

