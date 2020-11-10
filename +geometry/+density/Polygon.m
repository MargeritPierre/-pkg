classdef Polygon < pkg.geometry.density.Density
%POLYGON Density associated to a polygon
    
properties
    Points % polygon points
    Ce % polygon segments middle points
    Le % polygon segment lengths
end

methods
    function this = Polygon(points)
    % Class Constructor
        this.Points = points ;
        this.Ce = .5*(points(1:end,:) + points([2:end 1],:)) ;
        this.Le = sqrt(sum(diff(points([1:end 1],:),1,1).^2,2)) ;
    end
end

methods
    function h = evalAt(this,p)
    % Permute dimensions
        p = permute(p,[1 3 2]) ; % [nP 1 nCoord]
        ce = permute(this.Ce,[3 1 2]) ; % [1 nSeg nCoord]
        le = permute(this.Le,[2 1]) ; % [1 nSeg]
    % Edge distance
        dist = (sum((p-ce).^2,3)) ; % [nP nSeg]
        idist = 1./(dist+eps) ; % [nP nSeg] inverse of the distance
        sidist = sum(idist,2) ; % [nP 1] sum of inverse distances
    % Local density
        h = sum(le.*idist,2)./sidist ; % [nP 1]
    end
end
end

