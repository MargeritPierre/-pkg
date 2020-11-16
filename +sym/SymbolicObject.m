classdef (Abstract) SymbolicObject < handle ...
                                        & matlab.mixin.Heterogeneous ...
                                        & matlab.mixin.CustomDisplay 
%SYMBOLICOBJECT ABSTRACT base class for symbolic objects
    
properties
    Name char % name of the object
end

%% ABSTRACT PROPERTIES & METHODS
methods (Abstract , Access = protected)
    val = eval(this) ; % evaluate the function
end
methods (Abstract)
    fcn = diff(this,obj) % returns the function derivative w.r.t a symbolic object
end

%% CLASS CONSTRUCTOR

methods
    function this = SymbolicObject(name)
    % Class constructor
        if nargin==0 ; return ; end
        this.Name = name ;
    end
end

methods (Static , Access = protected)
    function obj = getDefaultScalarElement()
    % Return a default scalar symbolic object (for initialization purposes)
        obj = pkg.sym.Variable ;
    end
end

%% OPERATIONS
methods
    function fcn = uminus(this)
    % Return -fcn
        switch class(this)
            case 'pkg.sym.UMinus'
                fcn = this.Input ; % -(-f) = f
            case 'pkg.sym.Zero'
                fcn = this ; % -0 = 0
            otherwise
                fcn = pkg.sym.UMinus(this) ;
        end
    end
    
    function fcn = uplus(this)
    % Return -fcn
        fcn = this ;
    end
    
    function fcn = minus(this,obj)
    % Return -fcn
        switch class(obj)
            case 'pkg.sym.UMinus'
                fcn = pkg.sym.Plus(this,obj.Input) ; % A-(-B) = A+B
            otherwise
                fcn = pkg.sym.Minus(this,obj) ;
        end
    end
    
    function fcn = plus(this,obj)
    % Return -fcn
        switch class(obj)
            case 'pkg.sym.UMinus'
                fcn = pkg.sym.Minus(this,obj.Input) ; % A+(-B) = A-B
            otherwise
                fcn = pkg.sym.Plus(this,obj) ;
        end
    end
    
    function fcn = times(this,obj)
    % Return -fcn
        if isa(this,'pkg.sym.UMinus') && isa(obj,'pkg.sym.UMinus')
            fcn = pkg.sym.Times(this.Inputs,obj.Inputs) ; % (-A).*(-B) = A.*B
        elseif isa(this,'pkg.sym.UMinus')
            fcn = -pkg.sym.Times(this.Inputs,obj) ; % (-A).*B = -(A.*B)
        elseif isa(obj,'pkg.sym.UMinus')
            fcn = -pkg.sym.Times(this,obj.Inputs) ; % A.*(-B) = -(A.*B)
        else
            fcn = pkg.sym.Times(this,obj) ; % A.*B
        end
    end
    
    function fcn = mtimes(this,obj)
    % Return -fcn
        if isa(this,'pkg.sym.UMinus') && isa(obj,'pkg.sym.UMinus')
            fcn = pkg.sym.MTimes(this.Inputs,obj.Inputs) ; % (-A)*(-B) = A*B
        elseif isa(this,'pkg.sym.UMinus')
            fcn = -pkg.sym.MTimes(this.Inputs,obj) ; % (-A)*B = -(A*B)
        elseif isa(obj,'pkg.sym.UMinus')
            fcn = -pkg.sym.MTimes(this,obj.Inputs) ; % A*(-B) = -(A*B)
        else
            fcn = pkg.sym.MTimes(this,obj) ; % A*B
        end
    end
    
    function fcn = rdivide(this,obj)
    % Return -fcn
        fcn = pkg.sym.RDivide(this,obj) ;
    end
    
    function fcn = mrdivide(this,obj)
    % Return -fcn
        fcn = pkg.sym.MRDivide(this,obj) ;
    end
    
    function fcn = mldivide(this,obj)
    % Return -fcn
        fcn = pkg.sym.MLDivide(this,obj) ;
    end
    
    function fcn = power(this,obj)
    % Return -fcn
        fcn = pkg.sym.Power(this,obj) ;
    end
    
    function fcn = mpower(this,obj)
    % Return -fcn
        fcn = pkg.sym.MPower(this,obj) ;
    end
end

%% LOGICAL OPERATORS SEALED
methods (Sealed)
    function out = eq(A,B)
        out = eq@handle(A,B) ;
    end
    function out = ne(A,B)
        out = ne@handle(A,B) ;
    end
    function out = lt(A,B)
        out = lt@handle(A,B) ;
    end
    function out = le(A,B)
        out = le@handle(A,B) ;
    end
    function out = gt(A,B)
        out = gt@handle(A,B) ;
    end
    function out = ge(A,B)
        out = ge@handle(A,B) ;
    end
end

%% DISPLAY METHODS
methods
    function str = char(this)
    % Return a description of the symbolic object
        str = matlab.mixin.CustomDisplay.getClassNameForHeader(this) ;
    end
end
methods (Access = protected)
    function header = getHeader(this)
        header = ['symbolic ' char(this)] ;
        header = sprintf('\t%s\n',header) ;
    end
end

end

