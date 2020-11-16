classdef Zero < pkg.sym.Function

%% ABSTRACT PROPERTIES & METHODS
properties
    nInputs = NaN % number of inputs (can be NaN)
end

methods (Access = protected)
    function val = eval(this)
    % Evaluate the function at (numeric) places
        val = this.Inputs(1).eval() ;
        val = zeros(size(val),class(val)) ;
        for iii = 2:this.nInputs
            in = this.Inputs(iii).eval() ;
            val = val + zeros(size(in),class(in)) ;
        end
    end
end

methods
    function fcn = diff(this,obj)
    % Return the function derivative w.r.t a variable
        if nargin<2 ; obj = this.Variables(1) ; end
        fcn = pkg.sym.Zero(obj) ;
    end
end

end