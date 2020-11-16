classdef UMinus < pkg.sym.Operation

%% ABSTRACT PROPERTIES & METHODS
properties
    Symbols = {'-','',''}
    nInputs = 1 % number of inputs (can be NaN)
end

methods (Access = protected)
    function val = eval(this)
    % Evaluate the function at (numeric) places
        val = -this.Inputs.eval() ;
    end
end

methods
    function fcn = diff(this,obj)
    % Return the function derivative w.r.t a variable
        if nargin<2 ; obj = this.Variables(1) ; end
        fcn = -diff(this.Inputs,obj) ;
    end
end

end