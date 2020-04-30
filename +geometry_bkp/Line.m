classdef Line < pkg.geometry.Curve
    
    properties
        % [Start End]
        ControlPoints = NaN(2,3) ;
    end
    
    methods
        function this = Line(varargin)
        % Class construtor
            this = this@pkg.geometry.Curve() ;
            if nargin==1
                this.ControlPoints = varargin{1} ;
            elseif nargin==2
                this.ControlPoints = [varargin{1}(:)' ; varargin{2}(:)'] ;
            end
        end
    end
    
    methods
        function P = eval(this,t) 
        % Evaluate the curve at parameter t
            P = this.ControlPoints(1,:).*(1-t(:)) + this.ControlPoints(2,:).*t(:) ;
        end
    end
    
end