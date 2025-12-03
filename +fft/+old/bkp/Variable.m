classdef Variable < pkg.fft.old.Field
% VARIABLE containing gridded field data
    
    properties
        Data
    end
    
    methods
        function this = Variable(varargin)
        % Constructor
            this = this@pkg.fft.old.Field(varargin{:}) ;
        end
    end
end

