classdef EventData
% Contains all the data associated to the event data
    
    properties
        % TimeStamp
            TimeStamp = [] ;
        % Event
            Event = 'none'
        % Keyboard
            EventKey = ''
            ActiveKeys = {}
        % Mouse clicks
            EventMouseButton = ''
            ActiveMouseButtons = {}
        % Mouse scroll wheel
            ScrollAmount = NaN
        % Mouse position in figure [1 2] in pixels
            FigureMousePosition = NaN(1,2)
            FigureMouseMotion = NaN(1,2)
            FigureMouseMotionOrigin = NaN(1,2)
        % Mouse position in Axes [2 3] in data units
            AxesMousePosition = NaN(2,3)
            AxesMouseMotion = NaN(2,3)
            AxesMouseMotionOrigin = NaN(2,3)
        % Relative mouse position in Axes [1 2]
            RelMousePosition = NaN(1,2)
            RelMouseMotion = NaN(1,2)
            RelMouseMotionOrigin = NaN(1,2)
        % Objects bein under the mouse pointer
            MouseOn = gobjects(0)
    end
    
    methods
        function ED = EventData(VP,SRC,EVT)
        % Return all event data used by viewport interactions 
        % VP is a viewport
            % ARGUMENT PROCESSING
                if nargin==0 ; return ; end
            % INITIALIZE WITH THE PREVIOUS Viewport Event Data
                if ~isempty(VP.EventData) ; ED = VP.EventData(end) ; end
            % PROPERTIES ALWAYS CHANGING
                % Event
                    ED.TimeStamp = clock ;
                    ED.Event = EVT.EventName ;
                % Mouse Position
                    ED.AxesMouseMotion = VP.Axes.CurrentPoint-ED.AxesMousePosition ;
                    ED.AxesMousePosition = VP.Axes.CurrentPoint ;
                    ED.FigureMouseMotion = VP.Figure.CurrentPoint-ED.FigureMousePosition ;
                    ED.FigureMousePosition = VP.Figure.CurrentPoint ;
                    axPos = getpixelposition(VP.Axes) ;
                    rPos = (ED.FigureMousePosition-axPos(1:2))./axPos(3:4) ;
                    ED.RelMouseMotion = rPos-ED.RelMousePosition ;
                    ED.RelMousePosition = rPos ;
                    if ~strcmp(ED.Event,'WindowMouseMotion')
                        ED.AxesMouseMotionOrigin = ED.AxesMousePosition ;
                        ED.FigureMouseMotionOrigin = ED.FigureMousePosition ;
                        ED.RelMouseMotionOrigin = ED.RelMousePosition ;
                    end
                % Hitted object(s)
                    ED.MouseOn = hittest(VP.Figure) ; 
                    % Cull opaque objects (can happen on object deletion)
                        if strcmp(class(ED.MouseOn),'opaque') ; ED.MouseOn = VP.Axes ; end
                    % Cull non-viewport objects
                        if ~all(ismember(ED.MouseOn,[VP.Axes ; VP.Children(:)])) ;  ED.MouseOn = VP.Axes ; end
            % MOUSE EVENTS
                % Mouse clicks
                    ED.EventMouseButton = '' ;
                    if strcmp(ED.Event,'WindowMousePress')
                        ED.EventMouseButton = VP.Figure.SelectionType ;
                        ED.ActiveMouseButtons = flip(sort(unique([ED.ActiveMouseButtons {ED.EventMouseButton}]))) ;
                    elseif strcmp(ED.Event,'WindowMouseRelease')
                        ED.EventMouseButton = VP.Figure.SelectionType ;
                        if isempty(VP.Figure.SelectionType) || numel(ED.ActiveMouseButtons)==1
                            ED.ActiveMouseButtons = {} ;
                        else
                            ED.ActiveMouseButtons = setdiff(ED.ActiveMouseButtons,ED.EventMouseButton) ;
                        end
                    end
                % Scroll Wheel
                    if strcmp(ED.Event,'WindowScrollWheel')
                        ED.ScrollAmount = EVT.VerticalScrollCount ;
                    end
            % KEYBOARD EVENTS
                ED.EventKey = '' ;
                if strcmp(ED.Event,'WindowKeyPress')
                    ED.EventKey = EVT.Key ;
                    % Remove previous modifiers
                        ActiveKeys = setdiff(ED.ActiveKeys,EVT.Modifier,'stable') ;
                    % Add new keys
                        ActiveKeys = [ActiveKeys {EVT.Key}] ;
                    % Sort
                        ActiveKeys = sort(ActiveKeys) ;
                    % Add current modifiers
                        ActiveKeys = [EVT.Modifier ActiveKeys] ;
                    % Cull duplicates
                        ActiveKeys = unique(ActiveKeys,'stable') ;
                    if ~isequal(ED.ActiveKeys,ActiveKeys) % If the key combination is different
                        ED.ActiveKeys = ActiveKeys ;
                    else % The key is held
                        ED.Event = 'WindowKeyHold' ;
                    end
                elseif strcmp(ED.Event,'WindowKeyRelease')
                    ED.EventKey = EVT.Key ;
                    ED.ActiveKeys = setdiff(ED.ActiveKeys,ED.EventKey) ;
                end
        end
    end
    
end

