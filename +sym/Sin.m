classdef Sin < pkg.sym.Function

%% ABSTRACT PROPERTIES & METHODS
properties
    nInputs = 1 % number of inputs (can be NaN)
end

methods (Access = protected)
    function val = eval(this)
    % Evaluate the function at (numeric) places
        val = sin(this.Inputs.eval()) ;
    end
end

methods
    function fcn = diff(this,obj)
    % Return the function derivative w.r.t a variable
        if nargin<2 ; obj = this.Variables(1) ; end
        if obj==this % differentiation wrt itself=1
            fcn = pkg.sym.One(obj) ;
        else % chain rule
            fcn = pkg.sym.Cos(this.Inputs)*diff(this.Inputs,obj) ;
        end
    end
end

end