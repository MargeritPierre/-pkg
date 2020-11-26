classdef CATIA < pkg.ui.viewport.Interaction
% Pan, Rotate and Zoom the viewport with a CATIA V5 like behavior:
% everything is controlled by the middle and right mouse buttons
% middle btn: Pan
% middle then right btn: Rotate w.r.t axes center
% middle then right ON the right OFF : Zoom w.r.t axes center

    properties
        OrbitAmount = .25 % degrees/pixel
        RollAmount = 2 % ratio
        RollOrbitMargin = 0.1 % margin (relative to the axes size) in which roll mode is triggered
        ZoomAmount = 1.01 % zoomratio = zoomamount^motioninpixels
    end

    properties (AbortSet)
        CurrentMode char = 'none' % 'none', 'pan', 'orbit', 'roll' or 'zoom'
        PreviousFigurePointer = 'arrow' % viewport's figure pointer before mode change
        RollCircle = gobjects(0) ;
    end
    
%% CONSTRUCTOR/DESTRUCTOR
    methods
        function this = CATIA(varargin)
        % Class constructor
        end
        
        function delete(this)
        % Class destructor
        end
    end
    
%% INIT/UPDATE/REMOVE functions
    methods
        function init(this) 
        % Initialization function
        end
        
        function update(this) 
        % Updating function
            EventData = this.Viewport.EventData(end) ;
        % Is the event a compatible event ?
            if ~any(regexp(EventData.Event,'WindowMouse','once')) ... % it is not a mouse event
                || ~ismember('extend',EventData.ActiveMouseButtons) ... % the middle btn is not pressed
                || any(ismember({'normal','open'},EventData.ActiveMouseButtons)) % invalid btns are pressed
                this.CurrentMode = 'none' ;
                return ; 
            end
            switch EventData.Event
                case 'WindowMousePress'
                    switch EventData.EventMouseButton
                        case 'extend'
                            this.CurrentMode = 'pan' ;
                        case 'alt'
                            if any( EventData.RelMousePosition<this.RollOrbitMargin ...
                                    | EventData.RelMousePosition>1-this.RollOrbitMargin )
                                this.CurrentMode = 'roll' ;
                            else
                                this.CurrentMode = 'orbit' ;
                            end
                    end
                case 'WindowMouseRelease'
                    switch EventData.EventMouseButton
                        case 'extend'
                            this.CurrentMode = 'none' ;
                        case 'alt'
                            this.CurrentMode = 'zoom' ;
                    end
                case 'WindowMouseMotion'
                    switch this.CurrentMode
                        case 'pan'
                            this.Viewport.pan(EventData.AxesMouseMotionOrigin-EventData.AxesMousePosition) ;
                        case 'orbit'
                            orbitAngles = -2*pi/180*EventData.FigureMouseMotion*this.OrbitAmount ;
                            this.Viewport.rotate(orbitAngles) ;
                        case 'roll'
                            rollPos = EventData.RelMousePosition - [1;0].*EventData.RelMouseMotion ;
                            rollPolar = sum((rollPos-[0.5 0.5]).*[1 1i],2) ;
                            rollAngle = angle(rollPolar(2)/rollPolar(1))*this.RollAmount ;
                            this.Viewport.rotate([0 0 rollAngle]) ;
                        case 'zoom'
                            this.Viewport.zoom(this.ZoomAmount^EventData.FigureMouseMotion(2)) ;
                    end
            end
            %disp(['CATIA MODE: ',this.CurrentMode])
        end
        
        function remove(this) 
        % Removing function
            this.CurrentMode = 'none' ;
        end
    end
    
%% SET CURRENT MODE
methods
    function set.CurrentMode(this,mode)
    % Set the current interaction mode
        mode = lower(mode) ;
        if strcmp(this.CurrentMode,'none')
            this.PreviousFigurePointer = this.Viewport.Figure.Pointer ;
        end
        delete(this.RollCircle) ;
        switch mode
            case 'none'
                this.Viewport.Figure.Pointer = this.PreviousFigurePointer ;
            case 'pan'
                this.Viewport.Figure.Pointer = 'fleur' ;
            case 'orbit'
                this.Viewport.Figure.Pointer = 'circle' ;
            case 'roll'
                this.Viewport.Figure.Pointer = 'circle' ;
            case 'zoom'
                this.Viewport.Figure.Pointer = 'cross' ;
            otherwise
                error('unknown mode') ;
        end
        this.CurrentMode = mode ;
    end
    
    function initRollCircle(this)
    % Plot a circle when rolling mode is on
    end
end
    
%% VIEWPORT PROPERTIES CHANGE
    methods
        function onAxesChanged(this)
        % The viewport axes have changed
        end
    end
    
end

