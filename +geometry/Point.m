classdef Point < pkg.geometry.Geometry
% Geometrical object representing a point
    
    properties
        Position = [NaN NaN NaN]
    end
    
    methods
        function this = Point(varargin)
        % Class Constructor
            if nargin==1 % A [x y] or [x y z] vector is given
                this.Position = varargin{1}(:)' ;
            elseif nargin>1 % Coordinates are given separately
                this.Position = cat(2,varargin{:}) ;
            end
        end
    end
    
    methods
        function draw(this)
        % Drawing function
            draw@pkg.geometry.Geometry(this) ; % Empty the graphics list
            this.Patch.Vertices = this.Position ; 
        end
        
        function update(this)
        % Updating function
            if this.Visible
                this.Patch.Vertices = this.Position ;
            end
        end
    end
    
    methods
        function set.Position(this,pos)
        % Change the object's position
            this.Position = pos(:)' ;
            if numel(this.Position)<3 ; this.Position(3) = 0 ; end
            this.update ;
        end
    end
end

