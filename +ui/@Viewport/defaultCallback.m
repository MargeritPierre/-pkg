function defaultCallback(this,src,evt)

    % Infos
        % Camera Frame
            this.CameraFrame = this.getCameraFrame ;
        % Mouse Position
            this.MouseMotion = this.Axes.CurrentPoint-this.MousePosition ;
            this.MousePosition = this.Axes.CurrentPoint ;
        % Event
            EventName = evt.EventName ;
        % Hitted object
            this.MouseOn = hittest(this.Figure) ;
            if strcmp(class(this.MouseOn),'opaque') ; this.MouseOn = this.Axes ; end % Can happen on object deletion

    % On Mouse Event
        if strcmp(EventName,'WindowMousePress')
            this.ActiveMouseButton = this.Figure.SelectionType ;
            initSelectionRectangle(this) ;
        elseif strcmp(EventName,'WindowMouseRelease')
            this.ActiveMouseButton = 'none' ;
            this.SelectionRectangle.Visible = 'off' ;
            %this.Tooltip.Visible = 'off' ;
        end
        
    % On Scroll Wheel Event
        if strcmp(EventName,'WindowScrollWheel')
            ZoomAmount = 1.2^sign(evt.VerticalScrollCount) ;
            if any(isvalid(this.MouseOn)) && this.MouseOn(1)~=this.Axes && this.MouseOn(1)~=this.Figure
                center = mean([this.MouseOn(1).XData(:) this.MouseOn(1).YData(:) this.MouseOn(1).ZData(:)],1) ;
            else
                center = mean(this.MousePosition,1) ;
            end
            this.zoom(ZoomAmount,mean(this.MousePosition,1)) ;
        end
        
    % On Key Press/Release
        if strcmp(EventName,'WindowKeyPress')
            if ismember(evt.Key,evt.Modifier) % Only modifiers are pressed
                ActiveKey = strjoin(evt.Modifier,'+') ;
            else % Modifiers + Key
                ActiveKey = strjoin([evt.Modifier {evt.Key}],'+') ;
            end
            if ~strcmp(ActiveKey,this.ActiveKey) % If the key combination is different
                this.ActiveKey = ActiveKey ;
            else % The key is held
                %EventName = 'WindowKeyHold' ;
            end
        elseif strcmp(EventName,'WindowKeyRelease')
            this.ActiveKey = 'none' ;
        end
        
        
    % On GUI Resize
        if strcmp(EventName,'SizeChanged')
%             this.Axes.OuterPosition(3:4) = this.Figure.InnerPosition(3:4) ;
%             pixelsByUnit = mean(this.Axes.Position(3:4)./[range(this.Axes.XLim) range(this.Axes.YLim)]) ;
%             [Center,~] = this.getPosition() ;
%             Range = [this.Axes.Position(3) this.Axes.Position(4)]/pixelsByUnit ;
%             this.setPosition(Center,Range) ;
        end
        
    % Current Action Initialization
        if (   ... % IF Mouse or Keyboard action
               strcmp(EventName,'WindowKeyPress') ...
               || strcmp(EventName,'WindowKeyRelease') ...
               || strcmp(EventName,'WindowMousePress') ...
               || strcmp(EventName,'WindowMouseRelease') ...
            ) ...
            && ( ... % AND IF No key or no mouse btn active
               strcmp(this.ActiveMouseButton,'none') ...
               && strcmp(this.ActiveKey,'none') ...
            )
            this.CurrentAction = 'none' ;
            this.SelectionRectangle.Visible = 'off' ;
        end
        
    
    % Specific actions
        switch EventName
            case 'WindowMouseMotion'
                % Change the selection rectangle second corner (even if not
                % visible)
                    this.SelectionRectangle.Visible = 'off' ; 
                    this.updateSelectionRectangle ;
                switch this.CurrentAction
                    case 'Select'
                        this.SelectionRectangle.Visible = 'on' ; 
                        this.HighlightedObjects = this.objectsInsideSelectionRectangle ;
                        %this.SelectedObjects = inRect ;
%                     case 'AddToSelection'
%                         if ~strcmp(this.ActiveMouseButton,'none')
%                             this.SelectionRectangle.Visible = 'on' ;  
%                             inRect = objectsInsideSelectionRectangle ;
%                             [inRect.Selected] = deal(true) ;
%                         end
%                     case 'RemoveFromSelection'
%                         if ~strcmp(this.ActiveMouseButton,'none')
%                             this.SelectionRectangle.Visible = 'on' ;  
%                             inRect = objectsInsideSelectionRectangle ;
%                             [inRect.Selected] = deal(false) ;
%                         end
                    case 'Move'
%                         this.moveObjects(this.SelectedObjects,MouseMotion) ;
%                     case 'Rotate'
%                         % Compute the rotation angle
%                             dir = sign(angle(sum(MouseMotion.*[1 1i]))+eps) ;
%                             rot = dir*2*pi/180 ;
%                         % Apply
%                             this.rotateObjects(this.SelectedObjects,rot) ;
                    case 'Pan'
                        this.pan(this.MouseMotionOrigin-this.MousePosition) ;
                    case 'Rotate'
                        %this.rotate(this.MouseMotionOrigin-this.MousePosition) ;
                        this.rotate(-this.MouseMotion) ;
                    otherwise
                        this.HighlightedObjects = this.MouseOn ;
                end
            case 'WindowMousePress'
                this.MouseMotionOrigin = this.MousePosition ;
                switch this.Figure.SelectionType
                    case 'normal' % (left-click + no modifier)
                        switch this.MouseOn 
                            case this.Axes % User clicked somewhere in the canvas
                                this.CurrentAction = 'Select' ;
%                                 this.SelectedObjects = [] ;
                            otherwise % A Canvas object has been hitted
                                % Modify the selection if needed
%                                     if ~ismember(this.MouseOn,this.SelectedObjects) % A non-selected object is on the mouse pos.
%                                         this.SelectedObjects = this.MouseOn ;
%                                     end
                                % Initiate a motion
%                                     this.ActionOrigin = this.MousePosition ;
                                    switch this.ActiveKey
                                        case 'alt'
                                            this.CurrentAction = 'Rotate' ;
                                        otherwise
                                            this.CurrentAction = 'Move' ;
                                    end
                        end
                    case 'alt' % (right-click) or (control + left-click)
                        switch this.ActiveKey
                            case 'control' % (control + left-click)
%                                 switch this.MouseOn
%                                     case this % Click on the canvas, do nothing
%                                     otherwise % Click on an object, remove from selection
%                                         this.MouseOn.Selected = false ;
%                                 end
                            otherwise % (right-click) 
                                switch this.CurrentAction
                                    case 'Pan'
                                        this.CurrentAction = 'Rotate' ;
                                    otherwise
                                end
%                                 % Modify the selection if needed
%                                     if ~ismember(this.MouseOn,this.SelectedObjects) % A non-selected object is on the mouse pos.
%                                         this.SelectedObjects = this.MouseOn ; % returns empty if MouseOn is the Canvas
%                                     end
%                                 % Set the context menu
%                                     this.setContextMenu() ;
                        end
                    case 'open' % (double left-click)
%                             switch class(this.MouseHit)
%                                 case 'matlab.graphics.primitive.Patch'
%                                     if isa(this.MouseOn,'LINK.base.ui.canvasobjects.Node') ; this.addWire ; end
%                                     if isa(this.MouseOn,'LINK.base.ui.canvasobjects.Block') ; this.MouseOn.Peer.doubleClickCallback ; end
%                                 case 'matlab.graphics.primitive.Text'
%                                     this.MouseOn.Title.Editing = 'on' ; 
%                             end
                    case 'extend' % (center wheel) or (shift + left-click)
                        switch this.ActiveKey
                            case 'shift' % (shift + left-click)
%                                 switch this.MouseOn
%                                     case this % Click on the canvas, do nothing
%                                     otherwise % Click on an object, add to selection
%                                         this.MouseOn.Selected = true ;
%                                 end
                            otherwise % (scroll wheel click)
                                % Send object to workspace if possible
%                                     if this.MouseOn~=this
%                                         if ~isempty(this.MouseOn.Peer) && isvalid(this.MouseOn.Peer)
%                                             this.MouseOn.Peer.toWorkspace()
%                                         else
%                                             this.MouseOn.toWorkspace()
%                                         end
%                                         % Show object tooltip
%                                             this.setTooltip ;
%                                     end
                                % Initiate Panning
                                    this.CurrentAction = 'Pan' ;
%                                     this.ActionOrigin = MousePosition ;
                        end
                end
            case 'WindowKeyPress'
                switch this.ActiveKey
                    case 'control'
                        if ismember(this.CurrentAction,{'none','Select'})
                            this.CurrentAction = 'RemoveFromSelection' ;
                        end
                    case 'shift'
                        if ismember(this.CurrentAction,{'none','Select'})
                            this.CurrentAction = 'AddToSelection' ;
                        end
                    case 'escape'
                        switch this.CurrentAction
                            case 'none'
%                                 blocks = findType(this.Objects,'LINK.base.ui.canvasobjects.Block') ;
%                                 for comp = [blocks.Peer] ; comp.hideTab() ; end
                            otherwise
%                                 this.SelectedObjects = [] ;
                                this.SelectionRectangle.Visible = 'off' ;
                        end
                        this.CurrentAction = 'none' ;
                    case 'delete'
%                         delete(this.SelectedObjects) ;
                    case 'backspace'
%                         if ~isempty(this.SelectedObjects)
%                             set([this.SelectedObjects.Peer],'Active',false) ;
%                         end
                    case 'return'
%                         if ~isempty(this.SelectedObjects)
%                             set([this.SelectedObjects.Peer],'Active',true) ;
%                         end
                    case 'tab'
%                         if ~isempty(this.SelectedObjects)
%                             for obj = [this.SelectedObjects.Peer] ; obj.setPending() ; end
%                         end
                    case 'space'
%                         if numel(this.SelectedObjects)==1
%                             startUpdate(this.SelectedObjects.Peer) ;
%                         end
                    case 'control+a'
%                         this.SelectedObjects = this.Objects ;
                    case 'control+s'
%                         this.Session.save ;
                    case 'control+o'
%                         this.Session.load ;
                    case 'control+c'
%                         this.copyObjects(this.SelectedObjects) ;
                    case 'control+v'
%                         this.pasteObjects() ;
                    case 'control+d'
%                         this.dupplicateObjects(this.SelectedObjects) ;
                    case 'uparrow'
%                         if isempty(this.SelectedObjects)
%                             this.move([0 LINK.config.graphics.CanvasMotionIncrement]) ;
%                         else
%                             this.moveObjects(this.SelectedObjects,[0 LINK.config.graphics.CanvasMotionIncrement]) ;
%                         end
                    case 'downarrow'
%                         if isempty(this.SelectedObjects)
%                             this.move([0 -LINK.config.graphics.CanvasMotionIncrement]) ;
%                         else
%                             this.moveObjects(this.SelectedObjects,[0 -LINK.config.graphics.CanvasMotionIncrement]) ;
%                         end
                    case 'rightarrow'
%                         if isempty(this.SelectedObjects)
%                             this.move([LINK.config.graphics.CanvasMotionIncrement 0]) ;
%                         else
%                             this.moveObjects(this.SelectedObjects,[LINK.config.graphics.CanvasMotionIncrement 0]) ;
%                         end
                    case 'leftarrow'
%                         if isempty(this.SelectedObjects)
%                             this.move([-LINK.config.graphics.CanvasMotionIncrement 0]) ;
%                         else
%                             this.moveObjects(this.SelectedObjects,[-LINK.config.graphics.CanvasMotionIncrement 0]) ;
%                         end
                end
            case 'WindowKeyRelease'
                if ~strcmp(this.ActiveKey,'control') && strcmp(this.CurrentAction,'RemoveFromSelection')
                    this.CurrentAction = 'Select' ;
                end
                if ~strcmp(this.ActiveKey,'shift') && strcmp(this.CurrentAction,'AddToSelection')
                    this.CurrentAction = 'Select' ;
                end
        end
        
        
    % Mouse Pointer Appearance (let's have fun!!) <TODO>
        switch this.CurrentAction
            case 'Move'
                this.Figure.Pointer = 'arrow' ;
            case 'Pan'
                this.Figure.Pointer = 'fleur' ;
            case 'Rotate'
                this.Figure.Pointer = 'circle' ;
            otherwise
                this.Figure.Pointer = 'arrow' ;
        end
        
    % Display
        %if ~isempty(this.CurrentAction) && ~strcmp(this.CurrentAction,'none')
            %this.throw(this.CurrentAction) ; 
        %end
        if 0
            %disp(evt)
            %clc
            disp(newline)
            disp('Viewport Callback --------')
%             disp(['    EventName: ',EventName]) 
            disp(['    CurrentAction: ',this.CurrentAction]) 
%             disp(['    MousePos: ',mat2str(this.MousePosition)])
%             disp(['    MouseMotion: ',mat2str(this.MouseMotion)])
            disp(['    MouseMotionOrigin: ',mat2str(this.MouseMotionOrigin)])
%             disp(['    MouseOn: ',class(this.MouseOn)])
            disp(['    ActiveMouseButton: ',this.ActiveMouseButton])
%             disp(['    ActiveKey: ',this.ActiveKey])
%             disp(['    SelectedObjects: ',numel(this.SelectedObjects)])
%             disp(this.SelectedObjects)
%             disp(this.HighlightedObjects)
%             disp(this.SelectionRectangle)
        end
        
    % Draw
        %drawnow ;
        
    % Change Infos
        this.LastEvent = EventName ;








end
    
    
    
    