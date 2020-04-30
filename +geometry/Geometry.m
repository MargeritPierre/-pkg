classdef (Abstract) Geometry < handle & matlab.mixin.Copyable & matlab.mixin.Heterogeneous
% The base class of all geometrical objects
    
    properties % Geometries Hierarchy
        Parents = pkg.geometry.Geometry.empty
        Children = pkg.geometry.Geometry.empty
    end
    
    properties % User Interface
        % Graphical objects used for the reprezentation of the object
            Axes = gobjects(0) % <TODO> Multiple axes possible ?
            Patch = gobjects(0) ; % Graphical reprezentation
        % Object interaction
            Visible = false
            Selected = false
    end

    
%% CONSTRUCTOR/DESTRUCTOR
    methods
        function this = Geometry(varargin)
        % Class constructor
        end
        
        function delete(this)
        % Class destructor
            delete(this.Patch) ;
            delete(this.Children) ;
            this.Parents = [] ;
        end
    end
    
    
%% SET/GET PROPERTIES
    methods
        function set.Patch(this,graphics)
        % Change the list of graphics
            % Change the list 
                this.Patch = graphics ;
            % Apply the object's properties
                this.Axes = this.Axes ; 
                this.Visible = this.Visible ; 
        end
        
        function set.Parents(this,parents)
        % Change the object's parent(s)
            % Detach the object from its current parents
                for p = this.Parents(:)' ; p.Children(p.Children==this) = [] ; end
            % Attach the new parents
                this.Parents = parents ;
            % Add this as child to the new parents
                for p = this.Parents(:)' ; p.Children = unique([p.Children(:) ; this]) ; end
        end
        
        function set.Axes(this,ax)
        % Change the object's axes
            this.Axes = ax ;
            % Change object's graphic's parent
                set(this.Patch,'Parent',ax)
            % Change attached object's axes
                for c = this.Children(:)' ; c.Axes = ax ; end
        end
        
        function set.Visible(this,visible) 
        % Toggle the visibility of the object
            this.Visible = visible ;
            % Set the attached graphics visible or not
                switch this.Visible
                    case true
                        set(this.Patch,'Visible','on') ;
                    case false
                        set(this.Patch,'Visible','off') ;
                end
        end
        
        function set.Selected(this,selected)
        % Toggle the object selection
            this.Selected = selected ;
            this.setPatch ;
        end
    end
    
    
%% DRAW/UPDATE FUNCTIONS
    methods
        function draw(this)
        % Drawing function
            % Set the object's patch
                this.setPatch('init') ;
                this.setPatch('default') ;
            % Set default Axes
                if isempty(this.Axes) ; this.Axes = gca ; end
            % Make the object visible
                this.Visible = true ;
            % Draw attached objects
                for child = this.Children(:)' ; child.draw ; end
        end
        
        function update(this)
        % Updating function
            % Draw attached objects
                for child = this.Children(:)' ; child.update ; end
        end
    end
    
    
%% OBJECT APPEARANCE
    methods
        function setPatch(this,option)
        % Change the patch appearance
            % Choose the option if needed
                if nargin==1
                    if this.Selected
                        option = 'selected' ;
                    else
                        option = 'default' ; 
                    end
                end
            % Create the patch if needed
                if strcmp(option,'init')
                    delete(this.Patch) ;
                    this.Patch = patch(this.Axes,'Vertices',[NaN NaN NaN],'Faces',1) ;
                    this.Patch.UserData = this ;
                    this.Patch.PickableParts = 'visible' ;
%                     this.Patch.ButtonDownFcn = @(src,evt)this.buttonDownFcn ;
                end
            % If the patch exists
                if ~any(isvalid(this.Patch)) ; return ; end
            % Option-dependent properties
                switch option
                    case 'init'
                        this.Patch.Marker = 'o' ;
                        this.Patch.MarkerFaceColor = 'none' ;
                    case 'default'
                        this.Patch.LineWidth = 1 ;
                        this.Patch.MarkerSize = 5 ;
                        this.Patch.MarkerEdgeColor = 'k' ;
                    case 'selected'
                        this.Patch.LineWidth = 2 ;
                        this.Patch.MarkerSize = 7 ;
                        this.Patch.MarkerEdgeColor = 'r' ;
                    case 'highlighted'
                        this.Patch.LineWidth = 1.5 ;
                        this.Patch.MarkerSize = 6 ;
                end
        end
    end
    
    
%% OBJECT INTERACTIONS
    methods 
%         function buttonDownFcn(this)
%         % On a click on the object
%             if any(isvalid(this.Axes)) 
%                 f = this.Axes.Parent ;
%                 this.Axes.UserData.CurrentPoint = this.Axes.CurrentPoint ;
%                 f.WindowButtonUpFcn = @(src,evt)this.buttonUpFcn ; 
%                 f.WindowButtonMotionFcn = @(src,evt)this.mouseMotion(evt) ; 
%             end
%             disp('Down')
%             this.Selected = true ;
%         end
%         
%         function mouseMotion(this,evt)
%         % On mouse motion
%             %disp(evt)
%             motion = this.Axes.CurrentPoint-this.Axes.UserData.CurrentPoint ;
%             this.Axes.UserData.CurrentPoint = this.Axes.CurrentPoint ;
%             this.Position = this.Position + mean(motion,1) ;
%         end
%         
%         function buttonUpFcn(this)
%         % On a mouse button release
%             if any(isvalid(this.Axes))
%                 f = this.Axes.Parent ;
%                 f.WindowButtonUpFcn = [] ; 
%                 f.WindowButtonMotionFcn = [] ; 
%             end
%             disp('Up')
%             this.Selected = false ;
%         end
    end
    
end



