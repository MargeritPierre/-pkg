classdef Function < pkg.fft.old.Field
% FUNCTION a function acting on given gridded tensor
% fields/variables/functions

    properties
        OutputVariables pkg.fft.old.Variable
    end
    
    methods
        function this = Function(varargin)
        % Constructor
            this = this@pkg.fft.old.Field(varargin{:}) ;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
end
    
