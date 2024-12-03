function edgFcns = polylineEdgeFunctions(points,method)
%POLYLINEEDGEFUNCTIONS return the edge functions associated to a polyline
%defined by a series of points
% method: 'segments' or 'full'

    % Firts clean the list of points to remove too small segments
    tol = norm(range(points,1))*1e-6 ;
    dp = diff(points,1,1) ;
    tooSmall = sum(dp.^2,2)<tol^2 ;
    points = points([true ; ~tooSmall],:) ;
    
    
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
            edgFcns = {@(t)interp1(L/L(end),points,t(:),'linear','extrap')} ;
        otherwise % split edges with a given maximum angle
            edgvec = points(2:end,:) - points(1:end-1,:) ;
            edgvec = edgvec./sqrt(sum(edgvec.^2,2)) ;
            edgangle = acos(sum(edgvec(1:end-1,:).*edgvec(2:end,:),2)) ;
            splitidx = [1 ; find(edgangle>method)+1 ; size(points,1)] ;
            edgFcns = cell(1,numel(splitidx)-1) ;
            for ee = 1:numel(edgFcns)
                edgidx = splitidx(ee):splitidx(ee+1) ;
                edgFcns(ee) = pkg.geometry.levelset.polylineEdgeFunctions(points(edgidx,:),'full') ;
            end
    end
        
end

