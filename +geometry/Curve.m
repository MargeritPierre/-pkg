classdef (Abstract) Curve < pkg.geometry.Geometry
    
    properties (Abstract)
        % Control points defining the curve
        ControlPoints
    end
    
    methods (Abstract)
        P = eval(this,t) % Evaluate the curve at parameter t
    end
    
    methods
        function closed = isClosed(this,tol) 
        % Test if the curve is closed (the two end points coincides)
        end
        
        function planar = isPlanar(this,tol)
        % Test if the curve is planar
        end
        
        function L = getLength(this)
        % Compute approximately the length of the curve
        end
        
        function h = plot(this)
        % Plot the curve in the current axes
            ax = gca ;
            h = gobjects(0) ;
            P = this.eval(0:0.0001:1) ;
            h(end+1) = plot(ax,P(:,1),P(:,2),'k','linewidth',1.5) ;
                if size(P,2)>2 ; h(end).ZData = P(:,3) ; end
            h(end+1) = plot(ax,this.ControlPoints(:,1),this.ControlPoints(:,2),':k','linewidth',1) ;
                if size(this.ControlPoints,2)>2 ; h(end).ZData = this.ControlPoints(:,3) ; end
            h(end+1) = plot(ax,this.ControlPoints(:,1),this.ControlPoints(:,2),'.k','markersize',20) ;
                if size(this.ControlPoints,2)>2 ; h(end).ZData = this.ControlPoints(:,3) ; end
        end
    end
    
end