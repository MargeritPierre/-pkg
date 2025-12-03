classdef Operator < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
% OPERATOR on fields
    
    properties
        Input pkg.fft.Space
        Output pkg.fft.Space
    end
    
    methods
        function this = Operator(varargin)
        % Constructor
            if nargin==0 ; return ; end
            for arg = 1:2:nargin
                this.(varargin{arg}) = varargin{arg+1} ;
            end
        end

        function A = assemble(this)
        % Assemble the operator
        end

        function y = linearoperator(this,x)
        % return the linear operator corresponding to out = this.eval(in) ,
        % x = in.Data(:) and y = out.Data(:) ;
            in = pkg.fft.Field(this.Input,'Data',x) ;
            out = this.eval(in) ;
            y = out.Data(:) ;
        end
    end
        
    methods (Abstract)
        function out = eval(this,in)
        % Operator evaluation on pkg.fft.Fields
        end
    end
end

