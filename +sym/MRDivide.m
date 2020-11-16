classdef MRDivide < pkg.sym.Operation

%% ABSTRACT PROPERTIES & METHODS
properties
    Symbols = {'','/',''}
    nInputs = 2 % number of inputs (can be NaN)
end

methods (Access = protected)
    function val = eval(this)
    % Evaluate the function at (numeric) places
        val = this.Inputs(1).eval()/this.Inputs(2).eval() ;
    end
end

methods
    function fcn = diff(this,obj)
    % Return the function derivative w.r.t a variable
        if nargin<2 ; obj = this.Variables(1) ; end
        fcn = (diff(this.Inputs(1),obj)*this.Inputs(2) - this.Inputs(1)*diff(this.Inputs(2),obj))/this.Inputs(2).^2 ;
    end
end

end