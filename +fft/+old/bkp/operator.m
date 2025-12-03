classdef operator < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
%ABSTRACT OPERATOR on gridedd tensor fields (pkg.fft.old.field)
    
    properties
        Inputs pkg.fft.old.field
        Outputs pkg.fft.old.field
    end
    
    methods
        function this = operator(varargin)
        % Class constructor
            for arg = 1:2:numel(varargin)
                this.(varargin{arg}) = varargin{arg+1} ;
            end
        end

        function apply(this)
        % Function to be overrided by child classes
            if numel(this)>1 
                for oo = 1:numel(this) 
                    apply(this(oo)) ; 
                end
            end
            warning("Abstract operator function used !")
        end
    end
end

