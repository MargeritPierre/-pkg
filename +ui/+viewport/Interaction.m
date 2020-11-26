classdef Interaction < handle & matlab.mixin.Heterogeneous
%VIEWPORTINTERACTION base class for the viewport interaction classes
    
    properties
        % The parent viewport
        Viewport = pkg.ui.viewport.Viewport.empty
        % Active state
        Active = true
    end
    
%% CONSTRUCTOR/DESTRUCTOR
    methods
        function this = Interaction(varargin)
        % Class constructor
        end
        
        function delete(this)
        % Class destructor
            if ~isempty(this.Viewport) ; this.Viewport.removeInteraction(this) ; end
        end
    end
    
    methods (Static)
        function this = get(int)
        % Multiple object constructor (more flexible)
            this = pkg.ui.viewport.Interaction.empty ;
            if isempty(int) ; return ; end
            if isa(int,'pkg.ui.viewport.Interaction') ; this = int ; end
            switch class(int)
                case 'char' 
                    switch int
                        case 'all'
                            file = which('pkg.ui.viewport.Interaction') ;
                            [path,~,~] = fileparts(file) ;
                            folder = [path '\+interactions\*.m'] ;
                            classes = dir(folder) ;
                            classes = cellfun(@(c)c(1:end-2),{classes.name},'UniformOutput',false) ;
                            this = pkg.ui.viewport.Interaction.get(classes) ;
                        otherwise
                            interactionClass = ['pkg.ui.viewport.interactions.' int] ;
                            this = eval(interactionClass) ;
                    end
                case 'cell'
                    for ii = 1:numel(int) 
                        this(end+1) = pkg.ui.viewport.Interaction.get(int{ii}) ;
                    end
            end
            this = uniqueClass(this) ;
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
        
        function set.Active(this,state)
        % Toggle the interaction active state
            this.init ;
            this.Active = state ;
        end
    end
    
%% VIEWPORT PROPERTIES CHANGE
    methods
        function onAxesChanged(this)
        % The viewport axes have changed
        end
    end
    

%% OPERATIONS OVERRIDE
    methods (Sealed)
        function this = uniqueClass(this,varargin)
        % Keep only unique interactions
        % Avoid having two interactions of the same class in one viewport
            theseClasses = arrayfun(@(e)class(e),this,'UniformOutput',false) ;
            toRemove = [] ;
            for ii = 1:numel(this)
                int = this(ii) ;
                % Find other interactions with the same class
                    sameClass = find(ismember(theseClasses,class(int))) ;
                    sameClass = setdiff(sameClass,ii) ;
                % If needed, add to the remove list
                    if isempty(sameClass) 
                        continue ; 
                    else
                        toRemove(end+1) = ii ;
                        theseClasses{ii} = 'none' ;
                    end
            end
            this(toRemove) = [] ;
        end
        
        % Logical operators
        function out = eq(A,B) ; out = eq@handle(A,B) ; end
        function out = ne(A,B) ; out = ne@handle(A,B) ; end
        function out = lt(A,B) ; out = lt@handle(A,B) ; end
        function out = le(A,B) ; out = le@handle(A,B) ; end
        function out = gt(A,B) ; out = gt@handle(A,B) ; end
        function out = ge(A,B) ; out = ge@handle(A,B) ; end
    end
    
    
end

