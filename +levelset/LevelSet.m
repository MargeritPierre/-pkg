classdef LevelSet 
% The LevelSet class (signed distance function)
    

%% PROPERTIES
    properties
        % The function handle
        Function
        % The function domain bounds: [xmin ymin (zmin) ; xmax ymax (zmax)]
        BoundingBox
        % Edge fcn: a cell array of functon handles; P{i}(t) = EdgeFcns{i}(t), t in [0;1]
        EdgeFcns
        % Singular points
        Kinks
    end

    
%% CONSTRUCTOR / DESTRUCTOR
    methods
        function this = LevelSet(varargin)
        % Class Constructor
            if mod(nargin,2) ; error('wrong number of arguments') ; end
            for arg = 1:2:nargin-1
                this.(varargin{arg}) = varargin{arg+1} ;
            end
        end
    end
    
    
%% (STATIC) CONSTRUCTORS

    methods (Static)
        
        function ls = Segment(p1,p2)
        % LevelSet corresponding to a line
            ls = pkg.levelset.LevelSet() ;
            ls.Function = @(p)pkg.levelset.LevelSet.dsegment(p,p1,p2) ;
            ls.BoundingBox = [min(p1,p2) ; max(p1,p2)] ;
            ls.Kinks = [p1 ; p2] ;
            ls.EdgeFcns{1} = @(t)p1.*(1-t(:)) + p2.*t(:) ;
        end
        
        function ls = Polyline(points)
        % LevelSet corresponding to a polyline
            ls = pkg.levelset.LevelSet() ;
            ls.Function = @(p)pkg.levelset.LevelSet.dpolyline(p,points) ;
            ls.BoundingBox = [min(points,[],1) ; max(points,[],1)] ;
            ls.Kinks = points ;
            L = [0 ; cumsum(sqrt(sum(diff(points,1,1).^2,2)))] ;
            ls.EdgeFcns{1} = @(t)interp1(L/L(end),points,t,'linear','extrap') ;
        end
        
        function ls = Polygon(points)
        % LevelSet corresponding to a polygon
            ls = pkg.levelset.LevelSet() ;
            % Delete consecutive points that are too close
                dconsecutive = sqrt(sum((points-circshift(points,1,1)).^2,2)) ; 
                points(dconsecutive<eps,:) = [] ;
            ls.Function = @(p)pkg.levelset.LevelSet.dpolygon(p,points) ;
            ls.BoundingBox = [min(points,[],1) ; max(points,[],1)] ;
            ls.Kinks = points ;
            L = [0 ; cumsum(sqrt(sum(diff(points([1:end,1],:),1,1).^2,2)))] ;
            ls.EdgeFcns{1} = @(t)interp1(L/L(end),points([1:end,1],:),mod(t,1),'linear','extrap') ;
        end
        
        function ls = Rectangle(center,sides)
        % LevelSet corresponding to a rectangle
            ls = pkg.levelset.LevelSet() ;
            ls.Function = @(p)pkg.levelset.LevelSet.drectangle(p,center,sides) ;
            ls.Kinks = center + sides/2.*[-1 -1 ; 1 -1 ; 1 1 ; -1 1] ;
            ls.BoundingBox = ls.Kinks([1,3],:) ;
            points = ls.Kinks([1:end,1],:) ;
            L = [0 ; sides(1) ; sides(1)+sides(2) ; 2*sides(1)+sides(2) ; 2*sides(1)+2*sides(2)] ;
            ls.EdgeFcns{1} = @(t)interp1(L/L(end),points,mod(t,1),'linear','extrap') ;
        end
        
        function ls = Circle(center,radius)
        % LevelSet corresponding to a circle
            ls = pkg.levelset.LevelSet() ;
            ls.Function = @(p)pkg.levelset.LevelSet.dcircle(p,center,radius) ;
            ls.BoundingBox = center + radius*[-1 -1 ; 1 1] ;
            ls.Kinks = [] ;
            ls.EdgeFcns{1} = @(t)center + radius*[cos(2*pi*t(:)) sin(2*pi*t(:))] ;
        end
        
        function ls = Ellipse(center,semiaxes,rotation)
        % LevelSet corresponding to an ellipse
            ls = pkg.levelset.LevelSet() ;
            ls.Function = @(p)pkg.levelset.LevelSet.dellipse(p,center,semiaxes,rotation) ;
            ls.Kinks = [] ;
            ls.EdgeFcns{1} = @(t)center + [semiaxes(1)*cos(2*pi*t(:)) semiaxes(2)*sin(2*pi*t(:))]*[cos(rotation) sin(rotation) ; -sin(rotation) cos(rotation)] ;
            t = [-atan(semiaxes(2)/semiaxes(1)*tan(rotation)) ; atan(semiaxes(2)/semiaxes(1)/tan(rotation))] ;
            t = [t ; t+pi]/2/pi ;
            extrPts = ls.EdgeFcns{1}(t) ;
            ls.BoundingBox = [min(extrPts,[],1) ; max(extrPts,[],1)] ;
        end
    end
    
    
%% (STATIC) DISTANCE FUNCTIONS

    methods (Static)
        function d = dpoint(p,pos)
        % Distance to a point
            d = sqrt(sum((p-pos).^2,2)) ;
        end
        
        function [d,t] = dline(p,p1,p2)
        % Distance to an infinite line
            u = (p2-p1)./sqrt(sum((p2-p1).^2,2)) ;
            t = sum((p-p1).*u,2) ;
            d = sqrt(abs(sum((p-p1).^2,2) - t.^2)) ;
            if nargout>1 ; t = t./sqrt(sum((p2-p1).^2,2)) ; end
        end
        
        function d = dsegment(p,p1,p2)
        % Distance to a segment
            [d,t] = pkg.levelset.LevelSet.dline(p,p1,p2) ;
            dp1 = pkg.levelset.LevelSet.dpoint(p,p1) ;
            dp2 = pkg.levelset.LevelSet.dpoint(p,p2) ;
            d(t<=0) = dp1(t<=0) ;
            d(t>=1) = dp2(t>=1) ;
        end
        
        function d = dpolyline(p,points)
        % Distance to an open polyline
            points = permute(points,[3 2 1]) ;
            dsegments = pkg.levelset.LevelSet.dsegment(p,points(:,:,1:end-1),points(:,:,2:end)) ;
            d = min(dsegments(:,:),[],2) ;
        end
        
        function d = dpolygon(p,points)
        % (Signed) Distance to a polygon
            d = pkg.levelset.LevelSet.dpolyline(p,points([1:end,1],:)) ;
            d = d.*(0.5-inpolygon(p(:,1),p(:,2),points(:,1),points(:,2)))*2 ;
        end
        
        function d = drectangle(p,center,sides)
        % (Signed) Distance to a rectangle
            points = center + sides./2.*[-1 -1 ; 1 -1 ; 1 1 ; -1 1] ;
            d = pkg.levelset.LevelSet.dpolygon(p,points) ;
        end
        
        function d = dcircle(p,center,radius)
        % (Signed) Distance to a circle
            d = pkg.levelset.LevelSet.dpoint(p,center)-radius ;
        end
        
        function d = dellipse(p,center,semiaxes,rotation)
        % (Signed) Distance to an ellipse 
        % see https://wet-robots.ghost.io/simple-method-for-distance-to-ellipse/
            % We start by centering the ellipse on origin
                p = p-center ;
            % Then we rotate the frame of -rotation
                p = p*[cos(rotation) -sin(rotation) ; sin(rotation) cos(rotation)] ;
            % Any point pe on the ellipse is then :
                a = semiaxes(1) ; b = semiaxes(2) ;
                pe = @(t)[a*cos(t) b*sin(t)] ;
            % A point on the evolute is:
                pv = @(t)(b^2-a^2).*[-cos(t).^3./a sin(t).^3./b] ;
            % We push every point to the first quadrant
                p = abs(p) ;
            % And the initial guess is taken as follows:
                t = pi/4 ; atan2(a*p(:,2),b*p(:,1)) ;
            % Angle tolerance
                tol = 1e-6 ;
            % Maximum number of iterations
                maxIt = 5 ;
            % Newton descent
                dt = Inf ;
                it = 0 ;
                while it<maxIt && max(abs(dt))>tol
                    it = it+1 ;
                    %
                    x = pe(t) ;
                    e = pv(t) ;
                    %
                    r = x - e ;
                    q = p - e ;
                    %
                    R = sqrt(sum(r.^2,2)) ;
                    Q = sqrt(sum(q.^2,2)) ;
                    %
                    dc = R.*asin( (r(:,1).*q(:,2) - r(:,2).*q(:,1)) ./ (R.*Q) ) ;
                    dt = dc./sqrt(a^2 + b^2 - sum(x.^2,2)) ;
                    %
                    t = t + dt ;
                    t = min(pi/2,max(0,t)) ;
                end
            % Restablish to the four quadrants
                x = pe(t) ;
            % Absolute distance
                d = sqrt(sum((x-p).^2,2)) ;
            % Signed Distance
                d = d.*sign(p(:,1).^2/a^2+p(:,2).^2/b^2-1) ; 
        end
        
        function d = dunion(d1,d2)
        % Distance of union shapes
             d = min(d1,d2);
        end
        
        function d = ddiff(d1,d2)
        % Distance of substracted shapes
             d = max(d1,-d2);
        end
        
        function d = dintersect(d1,d2)
        % Distance of intersected shapes
             d = max(d1,d2);
        end
    end
    
    
%% BOOLEAN OPERATIONS
    methods
        function ls = or(ls1,ls2)
        % LevelSet function UNION
            ls = pkg.levelset.LevelSet() ;
            ls.Function = @(p)pkg.levelset.LevelSet.dunion(ls1.Function(p),ls2.Function(p)) ;
            ls.BoundingBox = [min(ls1.BoundingBox(1,:),ls2.BoundingBox(1,:)) ; max(ls1.BoundingBox(2,:),ls2.BoundingBox(2,:))] ;
            [ls.EdgeFcns,ls.Kinks] = mergeEdges(ls1,ls2) ;
            ls = ls.cleanContour ;
        end
        
        function ls = and(ls1,ls2)
        % LevelSet function INTERSECTION
            ls = pkg.levelset.LevelSet() ;
            ls.Function = @(p)pkg.levelset.LevelSet.dintersect(ls1.Function(p),ls2.Function(p)) ;
            ls.BoundingBox = [max(ls1.BoundingBox(1,:),ls2.BoundingBox(1,:)) ; min(ls1.BoundingBox(2,:),ls2.BoundingBox(2,:))] ;
            [ls.EdgeFcns,ls.Kinks] = mergeEdges(ls1,ls2) ;
            ls = ls.cleanContour ;
        end
        
        function ls = not(ls)
        % LevelSet function INVERSION
            ls.Function = @(p)-ls.Function(p) ; 
            ls.BoundingBox = Inf*[-1 -1 ; 1 1] ;
        end
        
        function ls = plus(ls1,ls2)
        % LevelSet function UNION
            ls = ls1 | ls2 ;
        end
        
        function ls = minus(ls1,ls2)
        % LevelSet function SUBSTRACTION
            ls = ls1 & ~ls2 ;
        end
    end
    
    
%% EGGES INTERSECTIONS

    methods
        function ti = intersectParam(this,edgfcn,t0)
        % Return the parameters t corresponding to the intersection of an
        % edge function with a levelset function
            if nargin<3 || isempty(t0) % Initialize the cross points
                t0 = linspace(0,1,1000)' ;
                p0 = edgfcn(t0) ;
                d0 = this.Function(p0) ;
                cross = sign(d0(1:end-1))~=sign(d0(2:end)) ;
                t0 = (t0([cross ; false]) + t0([false ; cross]))/2 ;
            end
            ti = t0 ;
            for ttt = 1:numel(t0)
                ti(ttt) = fzero(@(t)this.Function(edgfcn(t)),t0(ttt)) ;
            end
        end
        
        function [edgeFcns,kinks] = mergeEdges(ls1,ls2)
        % Process the edges of both levelsets to return the new contour
            edgeFcns = {} ;
            % Keep the initial singular points
                kinks = [ls1.Kinks ; ls2.Kinks] ;
            % For the two lvlst..
                LS = [ls1 ls2] ;
                for ll = 1:2
                    ls1 = LS(ll) ; ls2 = LS(3-ll) ;
                    % Cut each of the levelset's edges with the second's function
                        for ee = 1:numel(ls1.EdgeFcns)
                            % Find the intersections
                                ti = sort(ls2.intersectParam(ls1.EdgeFcns{ee})) ;
                            % Add the intersection points as singularities
                                kinks = [kinks ; ls1.EdgeFcns{ee}(ti)] ;
                            % Reparametrize edge functions
                                ti = unique([0 ; ti ; 1]) ;
                                for ii = 1:numel(ti)-1
                                    edgeFcns{end+1} = @(t)ls1.EdgeFcns{ee}((1-t(:))*ti(ii)+t(:)*ti(ii+1)) ;
                                end
                        end
                end
        end
        
        function ls = cleanContour(ls)
        % Clean the levelset edge functions and singular points
            % Delete invalid edges
                valid = true(size(ls.EdgeFcns)) ;
                for ee = 1:numel(ls.EdgeFcns)
                    [~,on] = ls.inside(ls.EdgeFcns{ee}((0:.1:1)')) ;
                    if any(~on) ; valid(ee) = false ; end
                end
                ls.EdgeFcns(~valid) = [] ;
            % Delete invalid kinks
                [~,on] = ls.inside(ls.Kinks) ;
                ls.Kinks(~on,:) = [] ;
                ls.Kinks = uniquetol(ls.Kinks,1e-6,'ByRows',1) ;
        end
    end


%% GEOMETRY UTILS
    methods
        function [in,on] = inside(this,P,dtol) 
        % Return a logical vector for points inside the domain / on the
        % domain boundary (up to dtol)
            if nargin<3 ; dtol = 1e-6*norm(range(this.BoundingBox,1)) ; end
            d = this.Function(P) ;
            in = d<dtol ;
            if nargin>1 ; on = in & d>-dtol ; end
        end
        
        function P = populate(this,density,distrib)
        % Fill the function domain with a distribution of points
        end
        
        function P = discretizeContour(this,dl)
        % Divide the contour edges by a target distance dl
            if nargin<2 ; dl = norm(range(this.BoundingBox,2))/100 ; end
            P = [] ;
            t0 = linspace(0,1,10)' ;
            for ee = 1:numel(this.EdgeFcns)
                P0 = this.EdgeFcns{ee}(t0) ;
                P = [P ; P0] ;
            end
        end
    end
    
    
%% VIZUALIZATION
    methods
        function h = plot(this,N)
        % Simple plot of the LevelSet
            % Function surface
                bbox = this.BoundingBox ;
                if isempty(bbox) ; bbox = [min(this.Kinks,[],1) ; min(this.Kinks,[],1)] ; end
                if isempty(bbox) ; bbox = [-1 -1 ; 1 1] ; end
                if norm(range(bbox,1))<eps ; bbox = bbox + [-1 -1 ; 1 1] ; end
                if nargin<2 ; N = 100*range(bbox,1)/norm(range(bbox,1)) ; end
                N = round(N).*[1 1] ;
                xx = linspace(bbox(1,1),bbox(2,1),N(1)) ;
                yy = linspace(bbox(1,2),bbox(2,2),N(2)) ;
                [XX,YY] = meshgrid(xx,yy) ;
                DD = reshape(this.Function([XX(:) YY(:)]),size(XX)) ;
                h = surf(XX,YY,DD,'facecolor','interp','edgecolor','none','FaceAlpha',0.5) ;
            % BoundingBox
                bboxP = [bbox(1,:) ; bbox(2,1) bbox(1,2) ; bbox(2,:) ; bbox(1,1) bbox(2,2)] ;
                h(end+1) = plot(bboxP([1:end,1],1),bboxP([1:end,1],2),':r') ;
            % Contour at d=0
                t = linspace(0,1,1000) ;
                for ee = 1:numel(this.EdgeFcns)
                    Pc = this.EdgeFcns{ee}(t) ;
                    h(end+1) = plot(Pc(:,1),Pc(:,2),'k') ;
                end
            % Kinks
                if ~isempty(this.Kinks)  
                    h(end+1) = plot(this.Kinks(:,1),this.Kinks(:,2),'.k','markersize',25) ;
                end
            % Display prop
                set(h,'Parent',gca) ; 
        end
    end



end

