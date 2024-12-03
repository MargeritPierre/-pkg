classdef LevelSet < matlab.mixin.Heterogeneous
% The LevelSet class (signed distance function)
    

%% PROPERTIES
    properties
        % The function handle
        Function
        % The function domain bounds: [xmin ymin (zmin) ; xmax ymax (zmax)]
        BoundingBox
        % Edge fcn: a cell array of functon handles; P{i}(t) = EdgeFcns{i}(t), t in [0;1]
        EdgeFcns
        % Surfacic mesh representation function (for quick intersection computation)
        SurfMeshFcn
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
    
    
%% (STATIC) DISTANCE FUNCTIONS

    methods (Static)
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
        
        function d = smoothmin(d1,d2,k,method)
        % Distance of union shapes with smooth minimums
        % see https://iquilezles.org/articles/smin/
            if nargin<3 ; k = 1 ; end
            if nargin<4 ; method = 'circular' ; end
            switch method
                case 'circular'
                    k = k./(1-sqrt(.5)) ;
                    h = max(k-abs(d1-d2),0)/k ;
                    d = min(d1,d2)-0.5*k*(1+h-sqrt(1-h.*(h-2))) ;
                otherwise
                    warning('Smooth minimum function not recognized. Applying min function.') ;
                    d = min(d1,d2) ;
            end
        end
    end
    
    
%% BOOLEAN OPERATIONS
    methods
        
        function ls = inversion(ls)
        % LevelSet function INVERSION
            ls.Function = @(p)-ls.Function(p) ; 
            ls.BoundingBox = Inf*([-1 ; 1]*ones(1,ls.nCoord)) ;
        end
        
        function ls = union(ls1,ls2)
        % LevelSet function UNION
            ls = pkg.geometry.levelset.LevelSet() ;
            ls.Function = @(p)pkg.geometry.levelset.LevelSet.dunion(ls1.Function(p),ls2.Function(p)) ;
            LS = [ls1 ls2] ; bbox = cat(1,LS.BoundingBox) ;
            ls.BoundingBox = [min(bbox,[],1) ; max(bbox,[],1)] ;
            [ls.EdgeFcns,ls.Kinks] = mergeEdges(ls1,ls2) ;
            ls.EdgeFcns = [ls.EdgeFcns ls1.intersectEdges(ls2)] ;
            ls = ls.cleanContour ;
        end
        
        function ls = substraction(ls1,ls2)
        % LevelSet function UNION (ls = ls1 & ~ls2) ;
            ls = pkg.geometry.levelset.LevelSet() ;
            ls.Function = @(p)pkg.geometry.levelset.LevelSet.ddiff(ls1.Function(p),ls2.Function(p)) ;
            ls.BoundingBox = ls1.BoundingBox ;
            [ls.EdgeFcns,ls.Kinks] = mergeEdges(ls1,ls2) ;
            ls.EdgeFcns = [ls.EdgeFcns ls1.intersectEdges(ls2)] ;
            ls = ls.cleanContour ;
        end
        
        function ls = intersection(ls1,ls2)
        % LevelSet function INTERSECTION
            ls = pkg.geometry.levelset.LevelSet() ;
            ls.Function = @(p)pkg.geometry.levelset.LevelSet.dintersect(ls1.Function(p),ls2.Function(p)) ;
            LS = [ls1 ls2] ; bbox = cat(3,LS.BoundingBox) ;
            ls.BoundingBox = [max(bbox(1,:,:),[],3) ; min(bbox(2,:,:),[],3)] ;
            [ls.EdgeFcns,ls.Kinks] = mergeEdges(ls1,ls2) ;
            ls.EdgeFcns = [ls.EdgeFcns ls1.intersectEdges(ls2)] ;
            ls = ls.cleanContour ;
        end
        
        function ls = merge(ls1,ls2,k)
        % Merge (union) two levelsets using smooth minimums
        % see https://iquilezles.org/articles/smin/
            ls = pkg.geometry.levelset.LevelSet() ;
            ls.Function = @(p)pkg.geometry.levelset.LevelSet.smoothmin(ls1.Function(p),ls2.Function(p),k) ;
            LS = [ls1 ls2] ; bbox = cat(1,LS.BoundingBox) ;
            ls.BoundingBox = [min(bbox,[],1) ; max(bbox,[],1)] + [-1;1]*k ;
            [ls.EdgeFcns,ls.Kinks] = mergeEdges(ls1,ls2) ;
            %ls.EdgeFcns = [ls.EdgeFcns ls1.intersectEdges(ls2)] ;
            ls = ls.cleanContour ;
        end
        
        % OVERRIDE BUILTIN FUNCTIONS
        function ls = or(ls1,ls2)
        % LevelSet function UNION
            ls = ls1.union(ls2) ;
        end
        
        function ls = and(ls1,ls2)
        % LevelSet function INTERSECTION
            ls = ls1.intersection(ls2) ;
        end
        
        function ls = not(ls)
        % LevelSet function INVERSION
            ls = ls.inversion() ;
        end
        
        function ls = plus(ls1,ls2)
        % LevelSet function UNION
            ls = ls1.union(ls2) ;
        end
        
        function ls = uplus(ls)
        % LevelSet function IDENTITY
        end
        
        function ls = minus(ls1,ls2)
        % LevelSet function SUBSTRACTION
            ls = ls1.substraction(ls2) ;
        end
        
        function ls = uminus(ls)
        % LevelSet function OPPOSITE
            ls = ls.inversion() ;
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
        
        function [edgeFcns,kinks] = intersectEdges(ls1,ls2)
        % Compute edge functions to be added when intersecting two (3D) levelsets
            edgeFcns = {} ; kinks = [] ;
            if ls1.nCoord<3 || ls2.nCoord<3 ; return ; end
        % One levelset will be discretized 
            % Base shape might have fastest meshBoundary overloads
%             isComplexShape = ismember({class(ls1) class(ls2)},{...
%                                         'pkg.geometry.levelset.LevelSet'...
%                                         }) ;
            isComplexShape = [isempty(ls1.SurfMeshFcn) isempty(ls2.SurfMeshFcn)] ;
            % Another option is the lvlst with the smallest discretization length
            [smallestLength,hasSmallestLength] = min([ls1.defaultDiscreteLength,ls2.defaultDiscreteLength]) ;
            % Choose the most adapted levelset to discretize
            LS = [ls1 ls2] ;
            if sum(isComplexShape)~=1 % if all ore none of the lvlst is a complex shape
                ls_d = LS(hasSmallestLength) ;
                ls_a = LS(3-hasSmallestLength) ;
            else
                ls_d = LS(~isComplexShape) ;
                ls_a = LS(isComplexShape) ;
            end
        % Discretize the chose levelset boundary
            bndmesh = ls_d.meshBoundary(smallestLength) ;
        % Compute its intersection curves with the other levelset
            cutmesh = bndmesh.cut(ls_a.Function).ON ;
            [~,bndCrv] = cutmesh.boundaryCurves ;
        % Create new edge functions with these curves
            maxAngle = pi/3 ; % edge split angle threshold
            for ee = 1:numel(bndCrv)
                polyfcn = pkg.geometry.levelset.polylineEdgeFunctions(cutmesh.Nodes(bndCrv{ee},:),maxAngle) ;
                for ff = 1:numel(polyfcn)
                    edgeFcns{end+1} = @(t)ls_d.toBoundary(polyfcn{ff}(t)) ; % re-project to boundary to reduce discretization error
                end
            end
        end
        
        function ls = cleanContour(ls,dtol)
        % Clean the levelset edge functions and singular points
            if nargin<3 ; dtol = 1e-3*norm(range(ls.BoundingBox,1)) ; end
            % Delete invalid edges
                valid = true(size(ls.EdgeFcns)) ;
                for ee = 1:numel(ls.EdgeFcns)
                    [~,on] = ls.inside(ls.EdgeFcns{ee}((0:.1:1)'),dtol) ;
                    if any(~on) ; valid(ee) = false ; end
                end
                ls.EdgeFcns(~valid) = [] ;
            % Delete invalid kinks
                if ~isempty(ls.Kinks)
                    [~,on] = ls.inside(ls.Kinks,dtol) ;
                    ls.Kinks(~on,:) = [] ;
                    ls.Kinks = uniquetol(ls.Kinks,1e-6,'ByRows',1) ;
                end
        end
    end
    
%% LEVELSET GRADIENT
methods
    function dgrad = gradient(this,p,normalize,geps)
    % Return the levelset gradient estimated at points p [nP nCoord]
    % normalize: boolean, normalize so that norm(grad)=1
    % geps: discrete step
    % grad: [nP nCoord]
        if nargin<3 ; normalize = true ; end
        if nargin<4 ; geps = norm(range(this.BoundingBox,1))*1e-6 ; end
    % Mean gradient over a quad/hex
        fd = @this.Function ; % distance function
        if isempty(p) ; dgrad = [] ; return ; end
        switch this.nCoord
            case 2
                sp = [-1 -1 ; 1 -1 ; 1 1 ; -1 1]/2 ; % [nG=4 nCoord] shifts
            case 3
                sp = [-1 -1 -1 ; 1 -1 -1 ; 1 1 -1 ; -1 1 -1 ; -1 -1 1 ; 1 -1 1 ; 1 1 1 ; -1 1 1]/2 ; % [nG=8 nCoord] shifts
        end
        pp = permute(p,[1 3 2]) + permute(sp,[3 1 2])*geps ; % [nP nG nCoord] all points
        d = reshape(fd(reshape(pp,[],size(p,2))),size(p,1),[]) ; % [nP nG] distance values
        dgrad = (d*sp)*(1/geps) ; % [nP nCoord] gradient
    % Normalize
        normGrad = sqrt(sum(dgrad.^2,2)) ;
        if normalize
            dgrad = dgrad./normGrad ;
        end
    % Check if the gradient is OK
        tooSmallNorm = normGrad<1e-1 ;
    % If not, shift the scheme by geps/2
        if any(tooSmallNorm) 
            warning('Vanishing gradient found !') ; 
            dgrad(tooSmallNorm,:) = this.gradient(p(tooSmallNorm,:)+geps/2,normalize,geps) ;
        end
    end
end

%% GEOMETRICAL OPERATIONS
    methods
        function this = move(this,v)
        % Move the levelset with a translation vector v
            this.Function = @(p)this.Function(p-v) ;
            this.BoundingBox = this.BoundingBox + v ;
            for ee = 1:numel(this.EdgeFcns)
                this.EdgeFcns{ee} = @(t)this.EdgeFcns{ee}(t) + v ;
            end
            this.Kinks = this.Kinks + v ;
        end
        
        function this = rotate(this,theta)
        % Rotate the levelset with an angle theta arounf the origin
            R = [cos(theta) sin(theta) ; -sin(theta) cos(theta)] ;
            this.Function = @(p)this.Function(p*R') ;
            this.BoundingBox = this.BoundingBox*R ;
            this.BoundingBox = [min(this.BoundingBox,[],1) ; max(this.BoundingBox,[],1)] ;
            for ee = 1:numel(this.EdgeFcns)
                this.EdgeFcns{ee} = @(t)this.EdgeFcns{ee}(t)*R ;
            end
            this.Kinks = this.Kinks*R ;
        end
        
        function this = offset(this,d)
        % Offset the levelset with a given distance
            this.Function = @(p)this.Function(p)-d ;
            this.BoundingBox = this.BoundingBox+[-1;1]*d ;
            this.EdgeFcns = {} ;
            this.Kinks = [] ;
        end
    end
%% GEOMETRY UTILS
    methods
        function n = nCoord(this) ; n = size(this.BoundingBox,2) ; end
        
        function [in,on] = inside(this,P,dtol) 
        % Return a logical vector for points inside the domain / on the
        % domain boundary (up to dtol)
            if nargin<3 ; dtol = defaultTolerance(this) ; end
        % first test in bbox
            bboxtol = this.BoundingBox + [-1;1]*dtol ;
            P = P(:,1:size(bboxtol,2),:,:,:,:,:,:) ;
            in = all( P>=bboxtol(1,:) & P<=bboxtol(2,:) ,2) ;
            if nargout>1 ; on = in ; end
        % Then use the distance function
            d = this.Function(P(in,:)) ;
            in(in) = in(in) & d<dtol ;
            if nargout>1 ; on(on) = on(on) & abs(d)<dtol ; end
        end
        
        function P = toBoundary(this,P,dtol)
        % Project a set of points to the closest levelset boundary
            if nargin<3 ; dtol = defaultTolerance(this) ; end
            maxIt = 100 ;
            outBND = true(size(P,1),1) ;
            it = 1 ;
            while it<maxIt && any(outBND)
                Pi = P(outBND,:) ;
                d = this.Function(Pi) ;
                g = this.gradient(Pi) ;
                P(outBND,:) = Pi - d.*g ;
                outBND(outBND) = abs(d)>dtol ; % convergence test
                it = it+1 ;
            end
        end
    end
    
%% GEOMETRY DISCRETIZATION
    methods
        function tol = defaultTolerance(this)
        % Default tolerance for geometrical operations
            tol = 1e-6*norm(range(this.BoundingBox,1)) ;
        end
        
        function h0 = defaultDiscreteLength(this)
        % Return a default discretization length
            h0 = max(range(this.BoundingBox,1))/30 ;
        end
        
        function P = populate(this,dx,distrib,fh,bnd)
        % Fill the function domain with a distribution of points
        % dx: spacing between points
        % distrib: point distribution
        %   - 'grid' : uniform square grid
        %   - 'iso' : isotropic grid of equilateral triangles
        %   - 'random' : random distribution
        % fh: (optional) local point spacing (edge length function)
        % bnd: populate levelset boundaries only
            L = range(this.BoundingBox,1) ; nCoord = numel(L) ;
            if nargin<2 ; dx = defaultDiscreteLength(this) ; end
            if nargin<3 ; distrib = 'grid' ; end
            if nargin<4 || isempty(fh) ; fh = @(P)dx ; end
            if nCoord~=2 && strcmp(distrib,'iso') 
                warning('No isotropic grid for dimension>2. Switching to regular grid') ;
                distrib = 'grid' ;
            end
            dx = dx(:)'.*ones(1,nCoord) ;
            % Initial distribution
                if strcmp(distrib,'iso') 
                    dx = [1 sqrt(3)/2].*dx ; 
                end
                switch distrib
                    case 'random'
                        nC = numel(L) ;
                        L = max(L) ;
                        N = ceil((L./dx(1))^nC) ;
                        P = rand(N,nC).*L + this.BoundingBox(1,:) ;
                    otherwise % 'grid' or 'iso'
                        % Initial grid
                            P = arrayfun(@colon ...
                                    ,this.BoundingBox(1,:)-dx ...
                                    ,dx ...
                                    ,this.BoundingBox(2,:)+dx ...
                                    ,'UniformOutput',false) ;
                            [P{:}] = ndgrid(P{:}) ;
                            P = cat(nCoord+1,P{:}) ;
                        % Shift if isotropic
                            if strcmp(distrib,'iso') 
                                P(:,1:2:end,1) = P(:,1:2:end,1) + dx(1)/2 ; 
                            end
                        % Reshape to a list of points
                            P = reshape(P,[],nCoord) ;
                end
            % Cull distribution if needed
                if nargin>3 
                    prob = prod(dx)./fh(P).^nCoord ; % rejection probability
                    pkeep = rand(size(P,1),1)<=prob ; % sampling
                    P = P(pkeep,:) ;
                end
            % Cull outside points
                P = P(this.inside(P),:) ;
            % If only boundaries needs to be populated
                if nargin>4 && bnd
                    on = abs(this.Function(P)) < .5*sqrt(nCoord)*fh(P) ;
                    P = P(on,:) ;
                    P = this.toBoundary(P) ;
                    P = uniquetol(P,this.defaultTolerance,'ByRows',true,'DataScale',1) ;
                end
        end
        
        function P = discretizeContour(varargin) % legacy compatibility
            P = discretizeEdges(varargin{:}) ; 
        end
        
        function P = discretizeEdges(this,dl,tol)
        % Divide the contour edges by a target distance dl
        % dl can be a scalar or a function handle dl = @(p)dl(P)
        % tol is used to cull dupplicate points: P = unique(P,tol*min(dl))
            if nargin<2 ; dl = defaultDiscreteLength(this) ; end
            if isa(dl,'function_handle') ; h0 = Inf ; else ; h0 = dl ; end
            if nargin<3 ; tol = 1/2 ; end
        % Initialize
            P = [] ;
            t0 = linspace(0,1,1000)' ;
        % For each edge
            for ee = 1:numel(this.EdgeFcns)
            % initial points
                P0 = this.EdgeFcns{ee}(t0) ;
            % segment lengths
                dL = sqrt(sum(diff(P0,1,1).^2,2)) ;
            % re-interpolation
                if isa(dl,'function_handle')
                    dl0 = dl(.5*(P0(1:end-1,:)+P0(2:end,:))) ; % local spacing
                    n = dL./dl0 ; % number of points in each interval
                    n = cumsum([0 ; n]) ; % cummulative number
                    n = (n./n(end)).*floor(n(end)) ;
                    if n(end)>0
                        t = interp1(n,t0,(0:n(end))') ; % parameters
                    else
                        t = [0;1] ;
                    end
                    h0 = min(h0,min(dl0(:))) ;
                else
                    L = [0 ; cumsum(dL,1)] ; % cummulative length
                    if L(end)>dl
                        Lt = linspace(0,L(end),floor(L(end)/dl)+1)' ;
                        t = interp1(L,t0,Lt,'linear') ;
                    else
                        t = [0;1] ;
                    end
                end
            % add to list
                P = [P ; this.EdgeFcns{ee}(t(:))] ;
            end
        % Add Kinks
            P = [P ; this.Kinks] ;
        % Cull duplicates
            P = uniquetol(P,tol*h0,'ByRows',true,'DataScale',1) ;
        end
        
        function skel = skeleton(this,varargin)
        % Return the levelset skeleton (2D only !)
            skel = pkg.geometry.Skeleton(this,varargin{:}) ;
        end
        
        function mesh = mesh(this,varargin)
        % Mesh the levelset
            mesh = pkg.geometry.mesh.distMesh(this,varargin{:}) ;
        end
        
        function mesh = meshBoundary(this,varargin)
        % Mesh the levelset boundary
            if nargin<3 && ~isempty(this.SurfMeshFcn) 
                if nargin==2 ; dx = varargin{1} ; 
                else ; dx = this.defaultDiscreteLength() ; 
                end
                mesh = this.SurfMeshFcn(dx) ;
            else
                mesh = this.mesh(varargin{:},'bnd_only',true) ;
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
                if nargin<2 
                    N = 1e5 ; L = range(bbox,1) ;
                    dx = (prod(L)/N)^(1/3) ;
                    N = round(L./dx) ; 
                end
                N = round(N).*ones(1,this.nCoord) ;
            % Discretization
                mesh = pkg.geometry.mesh.GridMesh(bbox,uint16(N)) ;
                X = mesh.Nodes ; D = this.Function(X) ;
                pl = plot(mesh,'VisibleFaces','all','VisibleEdges','none') ;
                switch this.nCoord
                    case 3 
                        pl.FaceAlpha = .05 ;
                    otherwise
                        pl.Deformation = [0 0 1].*D ;
                end
                pl.CData = D ; 
                pl.EdgeStyle = ':' ;
                h = pl.GraphicGroup ;
            % Contour at d=0
                t = linspace(0,1,1000) ;
                for ee = 1:numel(this.EdgeFcns)
                    Pc = this.EdgeFcns{ee}(t) ;
                    Pc(:,end+1:3) = 0 ; % 3D coordinates
                    h(end+1) = plot3(Pc(:,1),Pc(:,2),Pc(:,3),'k') ;
                end
            % Kinks
                if ~isempty(this.Kinks)  
                    Pk = this.Kinks ;
                    Pk(:,end+1:3) = 0 ; % 3D coordinates
                    h(end+1) = plot3(Pk(:,1),Pk(:,2),Pk(:,3),'.k','markersize',25) ;
                end
            % Display prop
                set(h,'Parent',gca) ; 
        end
    end



end

