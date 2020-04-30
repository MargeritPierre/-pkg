classdef Viewport < handle & matlab.mixin.Copyable
% This class defines a Viewport
% A viewport is a figure with interactions enabled
    
    properties % Graphical objects
        Figure = matlab.ui.Figure.empty
        Axes = matlab.graphics.axis.Axes.empty
    end
    
    properties % Interaction properties
        LastEvent = ''
        CurrentAction = 'none'
        ActiveKey = 'none'
        ActiveMouseButton = 'none'
        MousePosition = NaN(2,3)
        MouseMotion = NaN(2,3)
        MouseMotionOrigin = NaN(2,3)
        MouseOn = gobjects(0)
        CameraFrame = NaN(3,3)
    end
    
    properties % Objects
        SelectionRectangle = gobjects(0)
        BasePlanes = gobjects(0) 
    end
    
    properties (Dependent) % Lists
        HighlightedObjects
        SelectedObjects
    end


%% CONSTRUCTOR/DESTRUCTOR
    methods
        function this = Viewport(varargin)
        % Class Constructor
            if nargin==0
                this.Axes = gca ;
            elseif nargin==1
                this.Axes = varargin{1} ;
            end
            this.init ;
        end
        
        function delete(this)
        % Class Destructor
            % Reset properties to trigger eventual listeners
                this.reset ;
        end
    end
    
    
%% AXES/FIGURE/CAMERA
    methods
        function set.Axes(this,ax)
        % Change the viewport axes
            % Reset the current axes to default
                reset(this.Axes) ;
            % Set the new axes
                this.Axes = ax ;
            % Change the selection rectangle axes 
                set(this.SelectionRectangle,'Parent',ax) ;
            % Set the new axes behavior
                this.setAxesBehavior ;
            % Set the new figure (will transfer the callbacks)
                this.Figure = ax.Parent ;
        end
        
        function setAxesBehavior(this)
        % Set the axes properties
            reset(this.Axes) ;
            axis(this.Axes,'auto')
            axis(this.Axes,'normal')
            % User interactions
                this.Axes.Toolbar = [] ;
                this.Axes.Interactions = [] ;
            % Limits
                this.Axes.Clipping = 'off' ;
                this.Axes.XLimMode = 'manual' ;
                this.Axes.YLimMode = 'manual' ;
                this.Axes.ZLimMode = 'manual' ;
            % Grid
                grid(this.Axes,'off') ;
                set(this.Axes,'xtick',[],'ytick',[],'ztick',[]) ;
                set(this.Axes,'xcolor','none','ycolor','none','zcolor','none') ;
            % Appearance
                this.Axes.Projection = 'perspective' ;
                this.Figure.Renderer = 'painters' ;
                this.Axes.Box = 'on' ;
                this.Axes.BoxStyle = 'full' ;
            % Data aspect
                this.Axes.DataAspectRatio = [1 1 1] ;
            % Camera
                this.Axes.CameraViewAngle = 4 ;
        end
        
        function set.Figure(this,fig)
        % Change the object figure
            callbacks = this.getCallbacks ;
            % Apply the current callback state
                for cb = fieldnames(callbacks)'
                    fig.(cb{:}) = callbacks.(cb{:}) ;
                end
            % Delete all callbacks on the previous figure
                this.setCallbacks('all',[]) ;
            % Change
                this.Figure = fig ;
        end
        
        function frame = getCameraFrame(this)
        % Return the frame of the camera
            v3 = this.Axes.CameraPosition-this.Axes.CameraTarget ;
            v2 = this.Axes.CameraUpVector ;
            v1 = cross(v2,v3) ;
            v2 = cross(v3,v1) ;
            frame = [v1 ; v2 ; v3] ;
            frame = frame./sqrt(sum(frame.^2,2)) ;
        end
        
        function c = getAxesCenter(this)
        % Return the axes center point
            c = mean([this.Axes.XLim(:) this.Axes.YLim(:) this.Axes.YLim(:)],1) ;
        end
        
        function setAxesCenter(this,c)
        % Position the axes center (without changing the range)
            r = getAxesRange(this) ;
            this.Axes.XLim = c(1) + r(1)/2*[-1 1] ;
            this.Axes.YLim = c(2) + r(2)/2*[-1 1] ;
            this.Axes.ZLim = c(3) + r(3)/2*[-1 1] ;
        end
        
        function r = getAxesRange(this)
        % Return the axes center point
            r = range([this.Axes.XLim(:) this.Axes.YLim(:) this.Axes.YLim(:)],1) ;
        end
        
        function setAxesRange(this,r)
        % Position the axes center (without changing the range)
            c = getAxesCenter(this) ;
            this.Axes.XLim = c(1) + r(1)/2*[-1 1] ;
            this.Axes.YLim = c(2) + r(2)/2*[-1 1] ;
            this.Axes.ZLim = c(3) + r(3)/2*[-1 1] ;
        end
        
        function move(this,m)
        % Move the axes center
            this.setAxesCenter(this.getAxesCenter + m) ;
        end
        
        function moveTo(this,p)
        % Move the axes center
            this.setAxesCenter(p) ;
        end
        
        function pan(this,motion)
        % Pan the viewport with a motion (in the camera plane)
            motion = mean(motion,1) ;
            motion = motion*this.CameraFrame' ;
            motion = motion*(this.CameraFrame.*[1;1;0]) ;
            this.move(motion) ;
        end
        
        function rotate(this,motion)
        % Rotate the viewport with a motion (around the camera target)
            motion = mean(motion,1) ;
            motion = motion*this.CameraFrame' ;
            motion = motion./this.getAxesRange ;
            this.Axes.View = this.Axes.View + motion(1:2) * 150 ;
        end
        
        function zoom(this,factor,center)
        % Zoom the viewport by a factor
            if nargin<3 ; center = this.getAxesCenter ; end
            t = 1/factor ;
            this.setAxesCenter((1-t)*center + t*this.getAxesCenter) ;
            this.setAxesRange(t*this.getAxesRange) ;
        end
        
        function fitToObjects(this,objects)
        % Set the viewport so that all objects can be seen
            if nargin<2 ; objects = findobj(this.Axes,'-property','XData') ; end
            bbox = NaN(2,3) ;
            for obj = objects(:)'
                P = [obj.XData(:) obj.YData(:) obj.ZData(:)] ;
                if size(P)<3 ; P(:,3) = 0 ; end
                bbox = [min([bbox(1,:);P],[],1) ; max([bbox(2,:);P],[],1)] ;
            end
            if any(isnan(bbox(:))) || max(range(bbox,1))==0 ; bbox = [-1 -1 -1 ; 1 1 1] ; end
            this.setAxesRange(max(range(bbox,1))*[1 1 1]*2) ;
            this.setAxesCenter(mean(bbox,1)) ;
        end
    end
    
    
%% OBJECT SELECTION/HIGHLIGHTING
    methods
        function set.SelectedObjects(this,list)
        % Set the object selection
            % Remove all non-viewport childs from the list
                list(~ismember(list,this.Axes.Children)) = [] ;
            % Retrieve currently selected objects
                selected = this.SelectedObjects ;
            % Unselect all previously selected that are not in the list
                unselect = selected(~ismember(selected,list)) ;
                set(unselect,'Selected','off') ;
            % Select all objects in the list not already selected (UNHIGHLIGHT THEM BEFORE)
                list = list(~ismember(list,selected)) ;
                set(list,'SelectionHighlight','off') ;
                set(list,'Selected','on') ;
        end
        
        function selected = get.SelectedObjects(this)
        % Retrieve selected objects in the canvas
            selected = findobj(this.Axes,'Selected','on') ;
        end
        
        function set.HighlightedObjects(this,list)
        % Set the object in focus
            % Remove all non-viewport childs from the list
                list(~ismember(list,this.Axes.Children)) = [] ;
            % Retrieve currently highlighted objects
                highlighted = this.HighlightedObjects ;
            % Unhighlight all previously highlighted that are not highlighted
            % anymore
                unhighlight = highlighted(~ismember(highlighted,list)) ;
                set(unhighlight,'SelectionHighlight','off') ;
            % Highlight all newly highlighted objects THAT ARE NOT SELECTED
                list = list(~ismember(list,highlighted)) ;
                list = list(~ismember(list,this.SelectedObjects)) ;
                set(list,'SelectionHighlight','on') ;
        end
        
        function highlight = get.HighlightedObjects(this)
        % Retrieve selected objects in the canvas
            highlight = findobj(this.Axes,'SelectionHighlight','on') ;
        end
        
        function initSelectionRectangle(this)
        % On mouse button press
            this.updateSelectionRectangle('all') ;
        end
        
        function updateSelectionRectangle(this,points)
        % Update the selection rectangle position
            if nargin==1 ; points = 'second' ; end
            % Get the rectangle corners (in global coordinates)
                C2 = mean(this.MousePosition,1) ;
                switch points
                    case 'second'
                        C1 = [this.SelectionRectangle.XData(1) ...
                                this.SelectionRectangle.YData(1) ...
                                this.SelectionRectangle.ZData(1)] ;
                    case 'all'
                        C1 = C2 ;
                end
                Cg = [C1 ; C2] ;
            % Substract the camera Target
                Cg = Cg-this.Axes.CameraTarget ;
            % Project in the camera coordinates
                Cc = Cg*this.CameraFrame' ;
            % Compute the rectangle in the frame coordinate
                Rect = Cc(1,:).*[1 1 0 ; 0 1 0 ; 0 0 0 ; 1 0 0 ; 1 1 0] ...
                        + Cc(2,:).*[0 0 0 ; 1 0 0 ; 1 1 0 ; 0 1 0 ; 0 0 0] ;
            % Draw the rectangle in the [12] plane of the camera with an
            % origin at the camera target point (by default the axes
            % barycenter)
                Rect = this.Axes.CameraTarget + Rect*(this.CameraFrame.*[1;1;0]) ;
            % Update the display
                this.SelectionRectangle.XData = Rect(:,1) ;
                this.SelectionRectangle.YData = Rect(:,2) ;
                this.SelectionRectangle.ZData = Rect(:,3) ;
        end
        
        function inRect = objectsInsideSelectionRectangle(this)
        % Find the objects inside the rectangle
            inRect = gobjects(0) ;
            % Get the current bounding box in the camera frame
                R = [this.SelectionRectangle.XData([1;3])' this.SelectionRectangle.YData([1;3])' this.SelectionRectangle.ZData([1;3])'] ;
                R = R*this.CameraFrame' ;
                xmin = min(R(:,1)) ; xmax = max(R(:,1)) ;
                ymin = min(R(:,2)) ; ymax = max(R(:,2)) ;
            % Search for objects in this BBox
                objects = findobj(this.Axes.Children,'-property','XData') ;
                objects(objects==this.SelectionRectangle) = [] ;
                for obj = objects(:)'
                    P = [obj.XData(:) obj.YData(:) obj.ZData(:)] ;
                    if size(P,2)<3 ; P(:,3) = 0 ; end
                    P = P*this.CameraFrame' ;
                    isIn = P(:,1)>=xmin & P(:,1)<=xmax & P(:,2)>=ymin & P(:,2)<=ymax ; 
                    if any(isIn) ; inRect(end+1) = obj ; end
                end
        end
    end
    
    
%% INIT/UPDATE/RESET FUNCTIONS
    methods
        function init(this)
        % Initialization function
            this.reset ;
            this.setCallbacks ;
            this.initSelectionRectangle ;
        end
        
        function update(this)
        % Updating function
        end
        
        function reset(this)
        % Reset the viewport to an initial state
            % Objects
                delete(this.SelectionRectangle)
                this.SelectionRectangle = plot3(this.Axes,NaN,NaN,NaN,'-.k'...
                                                ,'linewidth',.5 ...
                                                ,'hittest','off' ...
                                                ,'PickableParts','none' ...
                                                ,'HandleVisibility','off' ...
                                                ) ;
                delete(this.BasePlanes) ;
                this.BasePlanes = patch(this.Axes ...
                                        ,'Vertices',[-1 -1 0 ; 1 -1 0 ; 1 1 0 ; -1 1 0 ... XY
                                                    ; 0 -1 -1 ; 0 1 -1 ; 0 1 1 ; 0 -1 1 ... YZ
                                                    ; -1 0 -1 ; 1 0 -1 ; 1 0 1 ; -1 0 1 ]/2 ... XZ
                                        ,'Faces',[1 2 3 4 ; 5 6 7 8 ; 9 10 11 12] ...
                                        ,'edgecolor','k' ...
                                        ,'facecolor','k' ...
                                        ,'facealpha',0.01 ...
                                        ,'linewidth',.5 ...
                                        ,'linestyle',':' ...
                                        ,'hittest','off' ...
                                        ,'PickableParts','none' ...
                                        ,'HandleVisibility','off' ...
                                        ) ;
                this.setAxesBehavior ;
                this.fitToObjects ;
            % Intercations
                this.LastEvent = '' ;
                this.CurrentAction = 'none' ;
                this.MouseOn = gobjects(0) ;
                this.MousePosition(:) = NaN ;
                this.MouseMotion(:) = NaN ;
                this.MouseMotionOrigin(:) = NaN ;
                this.ActiveMouseButton = 'none' ;
                this.ActiveKey = 'none' ;
                this.CameraFrame = NaN(3,3) ;
        end
    end
    
    
%% INTERACTION CALLBACKS
    methods
        function out = getCallbacks(this)
        % Return the callback state
            out = struct() ;
            if any(isvalid(this.Figure))
                for cb = {'SizeChangedFcn'...
                            ,'WindowKeyPressFcn'...
                            ,'WindowKeyReleaseFcn'...
                            ,'WindowButtonUpFcn'...
                            ,'WindowButtonDownFcn'...
                            ,'WindowButtonMotionFcn'...
                            ,'WindowScrollWheelFcn'}
                    out.(cb{:}) = this.Figure.(cb{:}) ;
                end
            end
        end
        
        function setCallbacks(this,actions,callback)
        % Set the callback state of the viewport
            if nargin<2 ; actions = 'all' ; end
            if nargin<3 ; callback = @(src,evt)this.defaultCallback(src,evt) ; end
            if any(isvalid(this.Figure))
                if ismember(lower(actions),{'all','size'})
                    this.Figure.SizeChangedFcn = callback ;
                end
                if ismember(lower(actions),{'all','keyboard'})
                    this.Figure.WindowKeyPressFcn = callback ;
                    this.Figure.WindowKeyReleaseFcn = callback ;
                end
                if ismember(lower(actions),{'all','mouse'})
                    this.Figure.WindowButtonUpFcn = callback ;
                    this.Figure.WindowButtonDownFcn = callback ;
                    this.Figure.WindowButtonMotionFcn = callback ;
                    this.Figure.WindowScrollWheelFcn = callback ;
                end
            end
        end
    end
    
    methods
        % The default Callback Behavior
        defaultCallback(this,src,evt)
    end
    
    
%% INTERACTIONS
    methods
    end


end