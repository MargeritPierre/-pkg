classdef Line < pkg.levelset.LevelSet
% Level set representing a Line

    properties
        Points
    end
    
    methods
        function this = Line(varargin)
        % Class construtor
            this = this@pkg.levelset.LevelSet() ;
            if nargin==2
                this.Points = [varargin{1}(:)' ; varargin{2}(:)'] ;
            elseif nargin==1
                this.Points = varargin{1}(:)' ;
            end
        end
    end
    
    methods
        function this = set.Points(this,pts)
        % Modify the points defining the line
            % Preliminary checks
                if size(pts,1)~=2 ; error('Wrong number of points: must be 2x2 or 2x3') ; end
                if norm(diff(pts,1))<eps ; error('The given points cannot be distinguished') ; end
            % Set the line points
                this.Points = real(pts) ;
                p1 = this.Points(1,:) ; p2 = this.Points(2,:) ;
            % Set the levelset properties
                this.Function = @(p)this.dline(p,p1,p2) ;
                this.EdgeFcns{1} = @(t)p1.*(1-t(:)) + p2.*t(:) ;
                this.BoundingBox = [] ;
                this.Kinks = [] ;
        end
    end
end

