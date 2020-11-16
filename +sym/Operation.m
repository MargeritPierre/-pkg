classdef (Abstract) Operation < pkg.sym.Function
% OPERATION a symbolic operation on symbolic objects
    
%% ABSTRACT PROPERTY
properties (Abstract)
    Symbols(1,3) cell % display : {presymbol symbol postsymbol}
end

%% COMMON PROPERTIES
methods
end

%% DISPLAY METHODS
methods
    function str = char(this)
    % Return a description of the symbolic function
        inputChar = arrayfun(@char,this.Inputs,'UniformOutput',false) ;
        isOp = arrayfun(@(i)isa(i,'pkg.sym.Operation'),this.Inputs) ;
        inputChar(isOp) = strcat('(',strcat(inputChar(isOp),')')) ;
        str = strjoin(inputChar,this.Symbols{2}) ;
        %str = ['(' str ')'] ;
        str = [this.Symbols{1} str this.Symbols{3}] ;
    end
end





end

