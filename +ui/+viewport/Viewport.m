classdef Viewport < handle & matlab.mixin.Copyable
% This class defines a Viewport
% A viewport is a figure with interactions enabled
    
    properties
        Axes matlab.graphics.axis.Axes = matlab.graphics.axis.Axes.empty
    end
    
    properties (SetAccess = private)
        Figure matlab.ui.Figure = matlab.ui.Figure.empty
    end
    
    properties (Dependent,SetAccess = protected)
        Children
    end
    
    properties % Viewport Interactions
        Interactions = pkg.ui.viewport.Interaction.empty ;
    end
    
    properties % Viewport Events
        EventBufferLength = 3 ;
        EventData pkg.ui.viewport.EventData = pkg.ui.viewport.EventData.empty
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
            this.Interactions(:) = [] ;
            this.setCallbacks('all','') ;
        end
    end
    
    
%% INTERACTIONS MANAGING
    methods
        function evalInteractions(this,fcn,interactions)
        % Make all (or given) viewport interactions execute a function
            if nargin<3 ; interactions = this.Interactions ; end
            for int = interactions(:)' 
                if int.Active ; fcn(int) ; end
            end
        end
        
        function set.Interactions(this,int)
        % Set the interactions of the viewport
            if isempty(int) ; int = pkg.ui.viewport.Interaction.empty ; end
            int = pkg.ui.viewport.Interaction.get(int) ;
            oldInt = this.Interactions ;
            % Set interactions
                this.Interactions = int(:) ;
            % Initialize new interactions
                newInt = int(~ismember(int,oldInt)) ;
                [newInt.Viewport] = deal(this) ;
                this.evalInteractions(@init,newInt) ;
            % Cull removed interactions
                remInt = oldInt(~ismember(oldInt,int)) ;
                this.evalInteractions(@remove,remInt) ;
        end
        
        function addInteraction(this,int)
        % Add an interaction to the viewport
            int = pkg.ui.viewport.Interaction.get(int) ;
            this.Interactions = [this.Interactions(:) ; int] ;
        end
        
        function removeInteraction(this,int)
        % Remove an interaction from the viewport
            this.Interactions(this.Interactions==int) = [] ;
        end
    end
    
    
%% INIT/UPDATE/RESET FUNCTIONS
    methods
        function init(this)
        % Initialization function
            this.EventData = pkg.ui.viewport.EventData.empty ;
            this.EventBufferLength = this.EventBufferLength ;
            this.setCallbacks ;
            this.evalInteractions(@init) ;
        end
        
        function update(this,src,evt)
        % Updating function
            % Make the viewport figure current (for keyboard capture)
                if gcf~=this.Figure ; return ; end
            % Update the event data
                newEvent = pkg.ui.viewport.EventData(this,src,evt) ;
                this.EventData = [this.EventData(2:end) ; newEvent] ;
            % Update the interactions
                t = tic ; 
                this.evalInteractions(@update) ; 
                %toc(t)
        end
        
        function reset(this)
        % Reset the viewport to an initial state
%             % Objects
%                 delete(this.SelectionRectangle)
%                 this.SelectionRectangle = plot3(this.Axes,NaN,NaN,NaN,'-.k'...
%                                                 ,'linewidth',.5 ...
%                                                 ,'hittest','off' ...
%                                                 ,'PickableParts','none' ...
%                                                 ,'HandleVisibility','off' ...
%                                                 ) ;
%                 delete(this.BasePlanes) ;
%                 this.BasePlanes = patch(this.Axes ...
%                                         ,'Vertices',[-1 -1 0 ; 1 -1 0 ; 1 1 0 ; -1 1 0 ... XY
%                                                     ; 0 -1 -1 ; 0 1 -1 ; 0 1 1 ; 0 -1 1 ... YZ
%                                                     ; -1 0 -1 ; 1 0 -1 ; 1 0 1 ; -1 0 1 ]/2 ... XZ
%                                         ,'Faces',[1 2 3 4 ; 5 6 7 8 ; 9 10 11 12] ...
%                                         ,'edgecolor','k' ...
%                                         ,'facecolor','k' ...
%                                         ,'facealpha',0.01 ...
%                                         ,'linewidth',.5 ...
%                                         ,'linestyle',':' ...
%                                         ,'hittest','off' ...
%                                         ,'PickableParts','none' ...
%                                         ,'HandleVisibility','off' ...
%                                         ) ;
%                 this.setAxesBehavior ;
%                 this.fitToObjects ;
        end
    end
    
    
%% AXES/FIGURE
    methods
        function set.Axes(this,ax)
        % Change the viewport axes
            % Reset the current axes to default
                reset(this.Axes) ;
            % Set the new axes
                this.Axes = ax ;
            % Set the new figure (will transfer the callbacks)
                if isempty(ax)
                    this.Figure = [] ;
                    return ;
                else
                    this.Figure = ax.Parent ;
                end
            % Set the new axes behavior
                this.setAxesBehavior ;
            % Notify interactions
                this.evalInteractions(@onAxesChanged) ;
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
                grid(this.Axes,'on') ;
                %set(this.Axes,'xtick',[],'ytick',[],'ztick',[]) ;
                %set(this.Axes,'xcolor','none','ycolor','none','zcolor','none') ;
            % Appearance
                %this.Axes.Projection = 'perspective' ;
                %this.Figure.Renderer = 'painters' ; % Vector graphics, but slow
                this.Axes.Box = 'on' ;
                this.Axes.BoxStyle = 'full' ;
            % Data aspect
                %this.Axes.DataAspectRatio = [1 1 1] ;
            % Camera
                %this.Axes.CameraViewAngle = 4 ;
        end
        
        function set.Figure(this,fig)
        % Change the object figure
        % Apply the current callback state to the new figure
            callbackSet = false ;
            if ~isempty(this.Figure) && ~isempty(fig)
                callbacks = this.getCallbacks ;
                for cb = fieldnames(callbacks)'
                    fig.(cb{:}) = callbacks.(cb{:}) ;
                end
                callbackSet = true ;
            end
        % Delete all callbacks on the previous figure
            if ~isempty(this.Figure)
                this.setCallbacks('all','') ;
            end
        % Change
            this.Figure = fig ;
        % Set default callback if needed
            if ~callbackSet && ~isempty(fig)
                this.setCallbacks() ;
            end
        end
    end
    
    
%% VIEWPORT CHILDREN
    methods
        function child = get.Children(this)
        % Return the viewport objects that are not part of interactions
            if ~isempty(this.Axes) && isvalid(this.Axes)
                child = this.Axes.Children ;
            else
                child = gobjects(0) ;
            end
            % <TODO> Cull invalid children
        end
    end
    
    
%% VIEWPORT MOTION CONTROL
    methods
        function frame = getCameraFrame(this)
        % Return the frame of the camera
            if isempty(this.Axes) || ~isvalid(this.Axes) ; frame = NaN(3,3) ; return ; end
            v3 = this.Axes.CameraPosition-this.Axes.CameraTarget ;
            v2 = this.Axes.CameraUpVector ;
            v1 = cross(v2,v3) ;
            v2 = cross(v3,v1) ;
            frame = [v1 ; v2 ; v3] ;
            frame = frame./sqrt(sum(frame.^2,2)) ;
        end
        
        function c = getAxesCenter(this)
        % Return the axes center point
            if isempty(this.Axes) || ~isvalid(this.Axes) ; c = NaN(1,3) ; return ; end
            c = mean([this.Axes.XLim(:) this.Axes.YLim(:) this.Axes.YLim(:)],1) ;
        end
        
        function r = getAxesRange(this)
        % Return the axes center point
            if isempty(this.Axes) || ~isvalid(this.Axes) ; r = NaN(1,3) ; return ; end
            r = range([this.Axes.XLim(:) this.Axes.YLim(:) this.Axes.YLim(:)],1) ;
        end
        
%         function setCameraFrame(this,fr)
%         % Set the frame of the camera while keeping the current target
%         % only the two first vecotrs of fr are used: the third is deduced
%             if isempty(this.Axes) || ~isvalid(this.Axes) ; return ; end
%             fr = fr./sqrt(sum(fr.^2,2)) ; % normalize
%             d = norm(this.Axes.CameraPosition-this.Axes.CameraTarget) ;
%             %this.Axes.CameraPosition = fr(3,:)
%             v3 = this.Axes.CameraPosition-this.Axes.CameraTarget ;
%             v2 = this.Axes.CameraUpVector ;
%             v1 = cross(v2,v3) ;
%             v2 = cross(v3,v1) ;
%             frame = [v1 ; v2 ; v3] ;
%             frame = frame./sqrt(sum(frame.^2,2)) ;
%         end
        
        function setAxesCenter(this,c)
        % Position the axes center (without changing the range)
            if isempty(this.Axes) || ~isvalid(this.Axes) ; return ; end
            r = getAxesRange(this) ;
            this.Axes.XLim = c(1) + r(1)/2*[-1 1] ;
            this.Axes.YLim = c(2) + r(2)/2*[-1 1] ;
            this.Axes.ZLim = c(3) + r(3)/2*[-1 1] ;
        end
        
        function setAxesRange(this,r)
        % Position the axes center (without changing the range)
            if isempty(this.Axes) || ~isvalid(this.Axes) ; return ; end
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
            CF = this.getCameraFrame ;
            motion = motion*CF' ;
            motion = motion*(CF.*[1;1;0]) ;
            this.move(motion) ;
        end

        function rotate(this,angles)
        % Rotate the viewport with angles (az,el,tilt)
        % angles are in radians
            if isempty(this.Axes) || ~isvalid(this.Axes) ; return ; end
            angles = angles/2/pi*180 ;
            camorbit(this.Axes,angles(1),angles(2)) ;
            if numel(angles)>2
                camroll(this.Axes,angles(3)) ;
            end
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
            if isempty(this.Axes) || ~isvalid(this.Axes) ; return ; end
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
    
    
%% VIEWPORTS EVENTS
    methods
        function set.EventBufferLength(this,len)
        % Change the event buffer length
            this.EventBufferLength = len ;
            if numel(this.EventData)>len
                this.EventData = this.EventData(end-len:end) ;
            else
                dummyEvt = pkg.ui.viewport.EventData() ;
                dummyEvts = repmat(dummyEvt,[len-numel(this.EventData) 1]) ;
                this.EventData = [dummyEvts(:) ; this.EventData(:)] ;
            end
        end
        
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
            if nargin<3 ; callback = @this.update ; end
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


end