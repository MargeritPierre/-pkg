classdef Variable < pkg.sym.SymbolicObject
% VARIABLE a symbolic variable

properties (Hidden)
    TempValue % value taken by the object temporarily
end

methods (Access = protected)
    function val = eval(this) 
    % Return the temporary value
        val = this.TempValue ;
    end
end

methods
    function fcn = diff(this,obj)
    % Variable differenctiation
        if nargin<2 ; obj = [] ; end
        if obj==this % differentiation wrt itself=1
            fcn = pkg.sym.One(obj) ;
        else % differentiation wrt other object=0
            fcn = pkg.sym.Zero(obj) ;
        end
    end
end

methods
    function str = char(this)
    % Return a description of the symbolic variable
        str = this.Name ;
    end
end

end

