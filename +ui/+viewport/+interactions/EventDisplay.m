classdef EventDisplay < pkg.ui.viewport.Interaction
% Display the event data on update

    properties
    % Declare here graphical objects, etc.
    end
    
%% CONSTRUCTOR/DESTRUCTOR
    methods
        function this = EventDisplay(varargin)
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
            disp(this.Viewport.EventData(end)) ;
        end
        
        function remove(this) 
        % Removing function
        end
    end
    
%% VIEWPORT PROPERTIES CHANGE
    methods
        function onAxesChanged(this)
        % The viewport axes have changed
        end
    end
    
end

