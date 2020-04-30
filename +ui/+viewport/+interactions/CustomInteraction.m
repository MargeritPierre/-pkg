classdef CustomInteraction < pkg.ui.viewport.Interaction
% Defines a custom interaction for a viewport
% To create a new interaction, find and replace 'CustomInteraction'

    properties
    % Declare here graphical objects, etc.
    end
    
%% CONSTRUCTOR/DESTRUCTOR
    methods
        function this = CustomInteraction(varargin)
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

