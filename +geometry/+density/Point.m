classdef Point < pkg.geometry.density.Density
%POINT Density associated to a Point
    
properties
    Center(1,:) % Position of the point
    Densities(:,2) % [distance density]
end

methods
    function this = Point(pt,dens)
    % Class Constructor
        this.Center = pt ;
        this.Densities = dens ;
    end
end

methods
    function this = set.Densities(this,dens)
    % Set function for densities
    % Get unique distance values (and sorted)
        [~,ia] = unique(dens(:,1)) ;
        dens = dens(ia,:) ;
    % Extrapolation close to the center
        if dens(1,1)>0 ; dens = [0 dens(1,2) ; dens] ; end
    % Set
        this.Densities = dens ;
    end
    
    function h = evalAt(this,p)
    % Evaluate the distance
        d = pkg.geometry.distance.point.toPoint(p,this.Center) ;
    % Interpolate the densities
        h = interp1(this.Densities(:,1),this.Densities(:,2),d,'linear',this.Densities(end,2)) ;
    end
end
end

