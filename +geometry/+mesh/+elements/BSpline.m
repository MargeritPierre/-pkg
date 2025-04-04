classdef BSpline < pkg.geometry.mesh.elements.AbstractElement
%BSpline Bi-Spline Elements

%% (MANDATORY) ELEMENT DEFINITION
% ELEMENT GEOMETRY
    properties (SetAccess = protected)
        % Node local Coordinates [nNodes nDims]
        NodeLocalCoordinates
        % The list of faces [nFaces nMaxNodesByFace]
        Faces = pkg.geometry.mesh.elements.ElementTable ;
        % The list of edges [nEdges nMaxNodesByEdge] 
        Edges = pkg.geometry.mesh.elements.ElementTable ;
    end
    
% SHAPE FUNCTIONS AND DERIVATIVES
    methods
        function C = eval_1D_At(this,u,q,p,der)
        % Evaluate the 1D BSpline functions using Cox-deBoor recursion formula
        % see https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
        % u is the local 1D coordinates in (0,1) [nE 1]
        % q is the scalar BSpline number of control points
        % p is the scalar BSpline polynomial degree
        % der is the (optionnal) derivation order (default 0)
        % C is the value of the BSplines coefficients at coordinates u [nE q]
            if nargin<5 ; der = 0 ; end
            if der>p ; C = zeros(numel(u),q) ; return ; end
        % Knot vector (assuming uniform spacing )
            ui = padarray((0:q-p)/(q-p),[0 p],'replicate','both') ;
        % Zero-th order basis (door) function
            C = u>=ui(1:end-1) & u<ui(2:end) ;
            C(u==1,end-p) = 1 ; % extremal values
        % Recursion formula
            for pp = 1:p-der
                kk = 1:size(C,2)-1 ; % knot spans
                C = (u-ui(kk))./max(ui(kk+pp)-ui(kk),eps).*C(:,1:end-1) ...
                    + (ui(kk+pp+1)-u)./max(ui(kk+pp+1)-ui(kk+1),eps).*C(:,2:end) ;
            end
        % Derivation
            if der>0
                for pp = p-der+1:p
                    kk = 1:size(C,2)-1 ; % knot spans
                    C = pp./max(ui(kk+pp)-ui(kk),eps).*C(:,1:end-1) ...
                        - pp./max(ui(kk+pp+1)-ui(kk+1),eps).*C(:,2:end) ;
                end
            end
        % Extrapolation by taylor expansion at order (p-der)
            outside = [u<0 u>1] ; 
            isout = any(outside,2) ;
            if ~any(isout) ; return ; end
            outside = outside(isout,:) ; % only keep meaningful info
            du = u(isout) - (outside*[0;1]) ; % deviation from the closest interval bound
            for o = 0:p-der
                C(isout,:) = C(isout,:) + (1/factorial(o))*(outside*this.eval_1D_At([0;1],q,p,der+o)).*(du.^o) ;
                %clf ; plot(u,C) ; drawnow
            end
        end
        
        % Evaluate the shape functions matrix N at local coordinates E
        % So that f(E) = N*f_nodes
        % E = [nE nDims] , N = [nE nNodes]
        function N = evalAt(this,E)
            N = ones(size(E,1),1) ;
            for d = 1:this.nDims
                C = this.eval_1D_At(E(:,d),this.Order(d),this.Degree(d)) ;
                N = repmat(N,1,size(C,2)).*repelem(C,1,size(N,2)) ;
            end
        end
        
        function DER = evalDerivativeAt(this,E,ORD,~)
        % Evaluate the partial derivatives of order ORD = [o1 o2 ...] (:[1 nDims]) of 
        % shape functions at given local coordinates E:[nRows nDims]
        % exemple: df(E)/(dx²dy) = DER(E,[2 1 0]) ;
            if nargin<3 ; ORD = [1 zeros(1,this.nDims-1)] ; end
            if numel(ORD)~=this.nDims ; error('Wrong derivation order argument (must be [1 nDims])') ; end
            if all(ORD==0) % No Derivative
                DER = this.evalAt(E) ;
            elseif any(ORD>this.Degree) % High-order derivatives vanish
                DER = zeros(size(E,1),this.nNodes) ;
            else % Valid derivatives
                DER = ones(size(E,1),1) ;
                for d = 1:this.nDims
                    C = this.eval_1D_At(E(:,d),this.Order(d),this.Degree(d),ORD(d)) ;
                    DER = repmat(DER,1,size(C,2)).*repelem(C,1,size(DER,2)) ;
                end
            end
        end
    end
    
% INTEGRATION QUADRATURE
    properties
        GaussIntegrationPoints % [nGaussIntPts nDims]
        GaussIntegrationWeights % [nGaussIntPts 1]
    end
    
%% CONSTRUCTOR / DESTRUCTOR
    properties
        Geometry = '1D' % element geometry: '1D', 'quad' or 'hex'
        Order = 5 % number of control points in each dimension [nDims 1]
        Degree = 3 % polynomial degree [nDims 1]
    end
    methods
        function this = BSpline(geo,ord,deg)
        % Constructor
        % Geometry
            if nargin>=1 ; this.Geometry = geo ; end
            switch this.Geometry
                case '1D' ; ndims = 1 ;
                case 'quad' ; ndims = 2 ;
                case 'hec' ; ndims = 3 ;
            end
        % Number of knots
            if nargin>=2 ; this.Order = ord ; end
            this.Order = this.Order(:).*ones(ndims,1) ;
        % Degree
            if nargin>=3 ; this.Degree = deg ; end
            this.Degree = this.Degree(:).*ones(ndims,1) ;
            if any(this.Degree>=this.Order)
                this.Degree = min(this.Degree,this.Order-1) ;
                warning('The BSpline degree has been adjusted in order to be lower than the corresponding order.') ;
            end
        % Node coordinates
            X = arrayfun(@colon,ones(ndims,1),ones(ndims,1).*this.Order,'UniformOutput',false) ;
            [X{:}] = ndgrid(X{:}) ;
            X = cat(ndims+1,X{:}) ;
            X = reshape(X-1,[],ndims)./(this.Order(:)'-1) ;
            this.NodeLocalCoordinates = X ;
        % Faces and edges
            IND = reshape(1:this.nNodes,[this.Order(:)' 1 1]) ;
            switch this.Geometry
                case '1D'
                    this.Edges = pkg.geometry.mesh.elements.ElementTable('Types',this,'NodeIdx',IND(:)') ;
                case 'quad'
                    this.Faces = pkg.geometry.mesh.elements.ElementTable('Types',this,'NodeIdx',IND(:)') ;
                    edgInd = padarray([1;2;1;2],[0 max(this.Order)],0,'post') ;
                        edgInd(1,2:this.Order(1)+1) = reshape(IND(:,1,1),1,[]) ;
                        edgInd(2,2:this.Order(2)+1) = reshape(IND(end,:,1),1,[]) ;
                        edgInd(3,2:this.Order(1)+1) = reshape(flip(IND(:,end,1)),1,[]) ;
                        edgInd(4,2:this.Order(2)+1) = reshape(flip(IND(1,:,1)),1,[]) ;
                    this.Edges = pkg.geometry.mesh.elements.ElementTable(...
                                        'Types', [ pkg.geometry.mesh.elements.BSpline('1D',this.Order(1),this.Degree(1)) ...
                                                   pkg.geometry.mesh.elements.BSpline('1D',this.Order(2),this.Degree(2)) ] ...
                                        ,'Indices',edgInd ) ;
                case 'hex'
            end
        end

        function delete(this)
        % Destructor
        end
    end

end

%% TEST FUNCTIONS
function test

%% CREATE A 1D BI-SPLINE ELEMENT AND PLOT ITS SHAPE FUNCTIONS & DERIVATIVES
ord = 6 ; deg = 3
elmt = pkg.geometry.mesh.elements.BSpline('1D',ord,deg) ;
clf ; plot(elmt,'shapefunctions')
%%
der = 1 ; 
%clf ; plot(elmt,'shapefunctions') ;

E = linspace(-0.5,1.5,100)' ; % extrapolation outside of [0 1]

N = elmt.evalAt(E) ;
clf ; plot(E,N) ;
if der<=0 ; return ; end

dN = elmt.evalDerivativeAt(E,der) ;
clf ; plot(E,dN) ;

ddN = N ;
for dd = 1:der
    ddN = [...
            diff(ddN([1 2],:),1,1) ; ...
            (ddN(3:end,:)-ddN(1:end-2,:))/2 ; ...
            diff(ddN(end-1:end,:),1,1) ...
          ]/mean(diff(E)) ;
end
set(gca,'ColorOrderIndex',1) ; plot(E,ddN,'.','markersize',15) ;


%% CREATE A 2D BI-SPLINE ELEMENT

elmt = pkg.geometry.mesh.elements.BSpline('quad',[1 1]*10,[4 4]) ;
der = [0 0] ;
%clf ; plot(elmt,'shapefunctions') ;

E = {linspace(0,1,100),linspace(0,1,100)} ;
[E{:}] = ndgrid(E{:}) ;
E = cat(numel(E)+1,E{:}) ;
N = reshape(elmt.evalDerivativeAt(reshape(E,[],size(E,ndims(E))),der),[size(E,1:ndims(E)-1) elmt.nNodes]) ;

cla ; axis square xy ; ax = gca ;
for nn = 27%:elmt.nNodes
    clr = mod(nn,size(ax.ColorOrder,1))+1 ;
    clr = ax.ColorOrder(clr,:) ;
    surf(E(:,:,1),E(:,:,2),N(:,:,nn),'facealpha',0.5,'facecolor',clr,'EdgeColor',clr) ;
end

%% 1D Node Localization
P = real(exp(2i*pi*(0:1/4:1-1e-9)').*[1 -1i]) ; n = size(P,1) ; p = min(n-1,2) ;
u = linspace(0,1,100)' ;
elmt = pkg.geometry.mesh.elements.BSpline('1D',n,p) ;
mesh = pkg.geometry.mesh.Mesh('Nodes',P,'Elems',pkg.geometry.mesh.elements.ElementTable('Types',elmt,'Indices',[1 1:elmt.nNodes])) ;
clf reset ; axis equal ; ax = gca ;
pl = patch('vertices',[elmt.evalAt(u)*P ; NaN NaN],'faces',1:numel(u)+1,'EdgeColor','b','FaceColor','none','LineWidth',2) ;
ma = patch('vertices',[NaN NaN],'faces',[1 1],'EdgeColor','b','FaceColor','b','Marker','.','MarkerSize',20) ;
poly = images.roi.Polyline('Parent',gca,'Position',P,'LineWidth',1) ;
addlistener(poly,'MovingROI',@(src,evt)set(pl,'vertices',[elmt.evalAt(u)*poly.Position ; NaN NaN])) ;
set(gcf,'WindowButtonMotionFcn',@(src,evt)set(ma,'vertices',elmt.evalAt(mesh.localize(ax.CurrentPoint(1,1:2),mesh.Elems,false,poly.Position))*poly.Position))

%% 1D BSPLINE MESH DERIVATIVE

elmt = pkg.geometry.mesh.elements.BSpline('1D',10,1) ;
mesh = pkg.geometry.mesh.Mesh('Nodes',elmt.NodeLocalCoordinates ...
                             ,'Elems',pkg.geometry.mesh.elements.ElementTable('Types',elmt,'NodeIdx',1:elmt.nNodes)) ;
clf ; plot(mesh,'VisibleNodes','all','NodeSize',25) ;

[E,ie] = mesh.localize(linspace(-.1,1.1,100)',mesh.Elems,true) ;
N = mesh.interpMat(E,ie) ;
D = mesh.diff2Mat(E,ie) ;
X = N*mesh.Nodes ;
F = X.^2 ;
Fn = N\F ;
dF = D{1}*Fn ;
clf ; plot(X,F) ; plot(X,N*Fn)
clf ; plot(X,dF)


%% INTEGRATION SCHEME ON BSPLINE MESHES

ord = 10 ; deg = 3 ;
elmt =  pkg.geometry.mesh.elements.BSpline('1D',ord,deg) ;

x = (0:ord-1)'*15/(ord-1) ;
mesh = elmt.mesh(x) ;

nQP = (ord-1)*deg ;

xiQP0 = .5/nQP+(0:nQP-1)'/nQP ;
xQP = (1-xiQP0).*x(1) + xiQP0.*x(end) ;
xiQP = mesh.localize(xQP) ;

elmt.GaussIntegrationPoints = xiQP ;
xiint = [0 ; .5*(xiQP(1:end-1)+xiQP(2:end)) ; 1] ;
Lxi = diff(xiint) ;
elmt.GaussIntegrationWeights = Lxi ;

range(x)

[ee,we,ie] = mesh.integration() ;
sum(we)
N = mesh.interpMat(ee,ie) ;

clf ;
plot(mesh,'VisibleNodes','all','NodeSize',15) ;
plot(N *mesh.Nodes,N*mesh.Nodes*0,'o')
plot(xQP,xQP*0,'+')


end
