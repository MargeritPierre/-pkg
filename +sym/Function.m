classdef (Abstract) Function < pkg.sym.SymbolicObject
% FUNCTION a symbolic function

%% ABSTRACT PROPERTIES & METHODS
properties (Abstract)
    nInputs % number of inputs (can be NaN)
end
    
%% COMMON PROPERTIES
properties
   Inputs pkg.sym.SymbolicObject % function "inputs"
   Variables pkg.sym.Variable % function "variables"
end

%% CLASS CONSTRUCTION
methods
    function this = Function(varargin)
    % Class constructor
        this.Name = this.defaultName ;
        if nargin==0 % no inputs, create default symbolic variables
            this.Inputs = defaultVariables(this) ;
        elseif isa(varargin{1},'pkg.sym.SymbolicObject') % symbolic variables given
            if ~isnan(this.nInputs) && nargin~=this.nInputs
                error('Wrong number of inputs given') ; 
            end
            this.Inputs = [varargin{:}] ;
            this.nInputs = numel(this.Inputs) ;
        else % invalid input
            error('Invalid input(s).') ;
        end
        this.Variables = this.findVariables ;
    end
    
    function name = defaultName(this)
    % Create the default name from the class's name
        name = regexprep(class(this),'pkg\.sym\.','') ;
    end
end

%% EVALUATION FUNCTION
methods
    function val = evalAt(this,varargin) 
    % Evaluate the function at non-symbolic values
    % Fill the variables's temporary value
        [this.Variables.TempValue] = deal(varargin{:}) ;
    % Execute the (abstract) evaluation function
        val = this.eval() ;
    % Clear all temporary values
        [this.Variables.TempValue] = deal([]) ;
    end
end

%% FUNCTIONS FOR INPUTS
methods
end

%% FUNCTION FOR VARIABLES
methods
    function n = nVars(this) 
    % Return the number of variables
        n = numel(this.Variables) ;
    end
    
    function vars = defaultVariables(this)
    % Set the default input number to 1
        if isnan(this.nInputs) ; this.nInputs = 1 ; end
    % Create default variables for the function
        vars(this.nInputs) = pkg.sym.Variable ;
    % default names: a,b,c,d,e,f,...
        defNames = arrayfun(@(n)char(n),96+(1:this.nInputs),'UniformOutput',false) ;
        [vars.Name] = deal(defNames{:}) ;
    end
    
    function vars = findVariables(this)
    % Find variables in the function inputs
    % Get the class of each input of the function
        inputClass = arrayfun(@class,this.Inputs,'UniformOutput',false) ;
    % Take variables apart
        isVar = strcmp(inputClass,'pkg.sym.Variable') ;
        vars = this.Inputs(isVar) ;
    % Recursively find variables in the function inputs
        if ~all(isVar)
            otherVars = arrayfun(@findVariables,this.Inputs(~isVar),'UniformOutput',false) ;
            vars = [cat(2,otherVars{:}) vars] ;
        end
    % Unique variables
        vars = unique(vars,'stable') ;
    end
end

%% DISPLAY METHODS
methods
    function str = char(this)
    % Return a description of the symbolic function
        className = matlab.mixin.CustomDisplay.getClassNameForHeader(this) ...
                    ... lower(this.Name) ...
                    ;
        inputChar = arrayfun(@char,this.Inputs,'UniformOutput',false) ;
        str = [className '(' strjoin(inputChar,',') ')'] ;
    end
end





end

