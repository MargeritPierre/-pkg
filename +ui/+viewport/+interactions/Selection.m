classdef Selection < pkg.ui.viewport.Interaction
% Allows selection of objects on the viewport

    properties
        CurrentAction = 'none' ;
        CurrentDomain  = 'point' ;
        MousePosition = NaN(2,3) ;
        MouseOn = gobjects(0) ;
    end
    
    properties
    % Declare here graphical objects, etc.
        SelectionRectangle = gobjects(0)
    end
    
    properties (Dependent) % Lists
        HighlightedObjects
        SelectedObjects
    end
    
    
%% CONSTRUCTOR/DESTRUCTOR
    methods
        function this = Selection(varargin)
        % Class constructor
        end
        
        function delete(this)
        % Class destructor
            delete(this.SelectionRectangle)
        end
    end
    
%% INIT/UPDATE/REMOVE functions
    methods
        function init(this) 
        % Initialization function
            this.CurrentAction = 'none' ;
            this.CurrentDomain = 'point' ;
            this.MousePosition(:) = NaN ;
            this.SelectedObjects = gobjects(0) ;
            this.HighlightedObjects = gobjects(0) ;
            delete(this.SelectionRectangle)
            if ~isempty(this.Viewport) && any(isvalid(this.Viewport)) ...
                    && ~isempty(this.Viewport.Axes) && isvalid(this.Viewport.Axes)
                this.SelectionRectangle = plot3(this.Viewport.Axes ...
                                                ,NaN*[1 1 1 1 1] ...
                                                ,NaN*[1 1 1 1 1] ...
                                                ,NaN*[1 1 1 1 1] ...
                                                ,'-.k'...
                                                ,'linewidth',.5 ...
                                                ,'hittest','off' ...
                                                ,'PickableParts','none' ...
                                                ,'HandleVisibility','off' ...
                                                ) ;
            end
            this.setSelection ;
        end
        
        function update(this) 
        % Updating function
            EventData = this.Viewport.EventData(end) ;
            this.MousePosition = EventData.AxesMousePosition ;
            this.MouseOn = EventData.MouseOn ;
            switch EventData.Event
                case 'WindowMousePress' 
                    % Begin selection on left-click
                        if isequal(EventData.ActiveMouseButtons,{'normal'}) 
                            this.initSelectionRectangle ;
                            switch EventData.MouseOn
                                case this.Viewport.Axes
                                    this.CurrentAction = 'init' ;
                                    this.CurrentDomain = 'rectangle' ;
                                otherwise
                                    this.CurrentAction = 'select' ;
                                    this.CurrentDomain = 'point' ;
                                    this.setSelection ;
                            end
                        end
                case 'WindowMouseMotion'
                    % Update highlighting on mouse motion
                        switch this.CurrentDomain
                            case 'point'
                                this.initSelectionRectangle ;
                                this.SelectionRectangle.Visible = 'off' ;
                            case 'rectangle'
                                switch this.CurrentAction
                                    case 'init'
                                        this.CurrentAction = 'highlight' ;
                                    case 'highlight'
                                        this.updateSelectionRectangle ;
                                        this.SelectionRectangle.Visible = 'on' ;
                                end
                        end
                        this.setSelection ;
                case 'WindowMouseRelease'
                    % Update selection on mouse release
                        switch this.CurrentAction
                            case 'init' % Reset the selection
                                this.CurrentAction = 'reset' ;
                                this.setSelection ;
                            case 'highlight'
                                if ismember(EventData.ActiveKeys,'shift')
                                    this.CurrentAction = 'addtoselection' ;
                                elseif ismember(EventData.ActiveKeys,'control')
                                    this.CurrentAction = 'removefromselection' ;
                                else
                                    this.CurrentAction = 'select' ;
                                end
                                this.setSelection ;
                        end
                        this.CurrentAction = 'highlight' ;
                        this.CurrentDomain = 'point' ;
                        this.SelectionRectangle.Visible = 'off' ;
            end
        end
        
        function remove(this) 
        % Removing function
            this.init ;
            delete(this.SelectionRectangle)
            this.Viewport.removeInteraction(this) ;
            this.Viewport = gobjects(0) ;
        end
    end
    
%% CURRENT ACTION SET
    methods
        function set.CurrentAction(this,action)
        % On the change of the current action
            this.CurrentAction = action ;
            %this.setSelection ;
        end
    end
    
%% VIEWPORT PROPERTIES CHANGE
    methods
        function onAxesChanged(this)
        % The viewport axes have changed
            set(this.SelectionRectangle,'Parent',this.Viewport.Axes) ;
        end
    end
    
%% SELECTION RECTANGLE DISPLAY
    methods
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
                CT = this.Viewport.Axes.CameraTarget ;
                Cg = Cg-CT ;
            % Project in the camera coordinates
                CF = this.Viewport.getCameraFrame ;
                Cc = Cg*CF' ;
            % Compute the rectangle in the frame coordinate
                Rect = Cc(1,:).*[1 1 0 ; 0 1 0 ; 0 0 0 ; 1 0 0 ; 1 1 0] ...
                        + Cc(2,:).*[0 0 0 ; 1 0 0 ; 1 1 0 ; 0 1 0 ; 0 0 0] ;
            % Draw the rectangle in the [12] plane of the camera with an
            % origin at the camera target point (by default the axes
            % barycenter)
                Rect = CT + Rect*(CF.*[1;1;0]) ;
            % Update the display
                this.SelectionRectangle.XData = Rect(:,1) ;
                this.SelectionRectangle.YData = Rect(:,2) ;
                this.SelectionRectangle.ZData = Rect(:,3) ;
        end
    end
    
    
%% OBJECT SELECTION/HIGHLIGHTING
    methods
        function set.SelectedObjects(this,list)
        % Set the object selection
            if isempty(this.Viewport) ; return ; end
            % Remove all non-viewport childs from the list
                list(~ismember(list,this.Viewport.Children)) = [] ;
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
            selected = findobj(this.Viewport.Children,'Selected','on') ;
        end
        
        function set.HighlightedObjects(this,list)
        % Set the object in focus
            if isempty(this.Viewport) ; return ; end
            % Remove all non-viewport childs from the list
                list(~ismember(list,this.Viewport.Children)) = [] ;
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
            highlight = findobj(this.Viewport.Children,'SelectionHighlight','on') ;
        end
        
        
        function setSelection(this,action,domain)
        % Set the selection data structure of objects
            if nargin<2 ; action = this.CurrentAction ; end
            if nargin<3 ; domain = this.CurrentDomain ; end
            % Depending on the current action
                %disp([action ' ' domain])
                switch action
                    case 'none'
                        apply = @(field,new,old)old.(field) ;
                    case 'reset'
                        apply = @(field,new,old)new.(field).*0 ;
                    case 'init'
                        apply = @(field,new,old)old.(field) ;
                    case 'highlight'
                        apply = @(field,new,old)max(new.(field)*1,(old.(field)==2)*2) ;
                    case 'select'
                        apply = @(field,new,old)new.(field)*2 ;
                    case 'addtoselection'
                        apply = @(field,new,old)max(new.(field)*2,old.(field)) ;
                    case 'removefromselection'
                        apply = @(field,new,old)min(~new.(field)*2,old.(field)) ;
                    otherwise
                        return ;
                end
            % Get the current bounding box in the camera frame
                % Global coordinates
                    if isvalid(this.SelectionRectangle)
                        BBox = [this.SelectionRectangle.XData([1;3])' this.SelectionRectangle.YData([1;3])' this.SelectionRectangle.ZData([1;3])'] ;
                    else
                        BBox = NaN(2,3) ;
                    end
                % Camera frame
                    CF = this.Viewport.getCameraFrame ;
                    BBox = BBox*CF' ;
                    BBox = sort(BBox(:,1:2),1) ;
                % Enhance bbox 
                    ran = this.Viewport.getAxesRange*CF' ;
                    BBox = BBox+[-1;1]*ran(1:2)/2*0.02 ;
            % Search for objects in this BBox
                for obj = this.Viewport.Children(:)'
                    % Get points of the object contained in the bbox
                        % Get 3D point coordinates
                            if isempty(obj.ZData) ; obj.ZData = obj.XData.*0 ; end
                            O = (obj.XData+obj.YData+obj.ZData)*0 ;
                            x = obj.XData+O ; y = obj.YData+O ; z = obj.ZData+O ;
                            P = [x(:) y(:) z(:)] ;
                            if size(P,2)<3 ; P(:,3) = 0 ; end 
                        % Project points in the camera frame
                            P = P*CF' ;
                        % Get members in the bbox
                            Pin = P(:,1)>=BBox(1,1) & P(:,1)<=BBox(2,1) & P(:,2)>=BBox(1,2) & P(:,2)<=BBox(2,2) ;
                    % Set selection data
                        newSel = struct() ;
                        newSel.Vertices = reshape(Pin,size(O)) ;
                        switch class(obj)
                            case 'matlab.graphics.chart.primitive.Surface'
                                newSel.Faces = newSel.Vertices(1:end-1,1:end-1) ...
                                                & newSel.Vertices(2:end,1:end-1) ...
                                                & newSel.Vertices(1:end-1,2:end) ...
                                                & newSel.Vertices(2:end,2:end) ;
                            case 'matlab.graphics.chart.primitive.Line'
                            case 'matlab.graphics.chart.primitive.Patch'
                        end 
                    % Retrieve previous selection data
                        if isempty(obj.UserData) || ~isstruct(obj.UserData)
                            obj.UserData = struct() ;
                        end
                        if isfield(obj.UserData,'Selection') && isstruct(obj.UserData.Selection)
                            oldSel = obj.UserData.Selection ;
                        else
                            oldSel = newSel ;
                        end
                    % Change the selection data
                        for fi = fieldnames(newSel)'
                            newSel.(fi{:}) = apply(fi{:},newSel,oldSel) ;
                        end
                        obj.UserData.Selection = newSel ;
                        %disp(newSel.Vertices)
                end
        end
    end
    
end

