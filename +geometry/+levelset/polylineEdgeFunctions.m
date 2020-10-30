function edgFcns = polylineEdgeFunctions(points,method)
%POLYLINEEDGEFUNCTIONS return the edge functions associated to a polyline
%defined by a series of points
% method: 'segments' or 'full'
    
    if nargin<2 ; method = 'segments' ; end

    switch method
        case 'segments' % One edge by segment of the polyline
            p1 = points(1:end-1,:) ; 
            p2 = points(2:end,:) ;
            nFcn = size(points,1)-1 ;
            edgFcns = cell(1,nFcn) ;
            for pp = 1:nFcn
                edgFcns{pp} = @(t)p1(pp,:).*(1-t(:)) + p2(pp,:).*t(:) ;
            end
        case 'full' % One edge reprezenting the full polyline
            L = [0 ; cumsum(sqrt(sum(diff(points,1,1).^2,2)))] ;
            edgFcns{1} = @(t)interp1(L/L(end),points,t(:),'linear','extrap') ;
    end
        
end

