function mesh = distMesh(lvlst,varargin)
%DISTMESH build a distance mesh from a pkg.geometry.levelset.LevelSet
% mesh = distMesh(lvlst) automatic meshing
% mesh = distMesh(lvlst,h0) uniform edge length h0 (scalar)
% mesh = distMesh(lvlst,'Name','Value',...) custom parameters

% Process input
    if ~isa(lvlst,'pkg.geometry.levelset.LevelSet')
        error('The first argument MUST be a pkg.geometry.levelset.LevelSet object') ;
    end
    if nargin==2 % given edge length
        varargin = {'h0',varargin{1}} ;
    end
    args = struct(varargin{:}) ; 

% Default input arguments 
    % Distance function
        if isfield(args,'fd') ; fd = args.fd ;
        else ; fd = @lvlst.Function ;
        end
    % MESH DENSITY
        % Initial (minimum) edge length 
            if isfield(args,'h0') ; h0 = args.h0 ;
            else ; h0 = lvlst.defaultDiscreteLength() ;
            end
        % Space-dependent edge length 
            if isfield(args,'fh') ; fh = args.fh ;
            else ; fh = @(p)h0*(ones(size(p,1),1)) ;
            end
        % Initial node distribution
            if isfield(args,'p0') ; p0 = args.p0 ;
            else ; p0 = lvlst.populate(h0,'iso',fh) ;
            end
        % Fixed nodes
            if isfield(args,'pfix') ; pfix = args.pfix ;
            else ; pfix =  lvlst.discretizeContour(fh) ... % variable-spaced contour points
                        ... lvlst.discretizeContour(h0) ... % uniform spaced contour points
                        ... lvlst.Kinks ... % only kink points
                        ... [] ... % no points
                        ;
            end
    % RELAXATION PARAMETERS
        % Max number of iterations
            if isfield(args,'maxCount') ; maxCount = args.maxCount ;
            else ; maxCount = 100 ;
            end
        % Node displacement tolerance (convergence criterion, relative to local edge length)
            if isfield(args,'dptol') ; dptol = args.dptol;
            else ; dptol = 0.01 ;
            end
        % Re-triangulation tolerance (relative to local edge length)
            if isfield(args,'ttol') ; ttol = args.ttol ;
            else ; ttol = 0.01 ;
            end
        % Minimum triangle quality
            if isfield(args,'qmin') ; qmin = args.qmin ;
            else ; qmin = 0.1 ; 
            end
        % Spring L0 relative to bar length
            if isfield(args,'Fscale') ; Fscale = args.Fscale ;
            else ; Fscale = 1.2 ; 
            end
        % Implicit time increment (relaxation damping)
            if isfield(args,'deltat') ; deltat = args.deltat ;
            else ; deltat = 1 ;
            end
        % Maximum distance function allowed for Nodes (relative to local edge length)
            if isfield(args,'p_dmax') ; p_dmax = args.p_dmax ;
            else ; p_dmax =  0.5 ; 0.001 ;
            end
        % Maximum distance function allowed for Triangle centroids (relative to local edge length)
            if isfield(args,'t_dmax') ; t_dmax = args.t_dmax ;
            else ; t_dmax = -0.2 ;
            end
        % Too SHORT edges threshold (rel. to local edge length)
            if isfield(args,'tooShortThrs') ; tooShortThrs = args.tooShortThrs ;
            else ; tooShortThrs = 0.7 ; 0.49 ;
            end
        % Too LONG edges threshold (rel. to local edge length)
            if isfield(args,'tooLongThrs') ; tooLongThrs = args.tooLongThrs ;
            else ; tooLongThrs = 1.3 ; 2.01 ;
            end
        % Constraint boundary nodes to be on the levelset edge
            if isfield(args,'bndCons') ; bndCons = args.bndCons ;
            else ; bndCons = true ;
            end
        % Discrete gradient estim. step
            if isfield(args,'geps') ; geps = args.geps ;
            else ; geps = 1e-2*h0 ; 
            end
    % MESH PREVIEW
        % Plot on every mesh change ?
            if isfield(args,'debug') ; debug = args.debug ;
            else ; debug = false ; 
            end
        % Mesh relaxation preview
            if isfield(args,'showMesh') ; showMesh = args.showMesh ;
            else ; showMesh = false || debug ; 
            end
        % Update plot frequency
            if isfield(args,'plotFreq') ; plotFreq = args.plotFreq ;
            else ; plotFreq = Inf ; 
            end
    % DISPLAY INFOS
        % Display info frequency
            if isfield(args,'displayFreq') ; displayFreq = args.displayFreq ;
            else ; displayFreq = 0 ; 
            end

% Remove points outside the geometry
    p = p0(fd(p0)<p_dmax.*fh(p0),:) ;
    
% Append fixed nodes
    if ~isempty(pfix), p=setdiff(p,pfix,'rows'); end     % Remove duplicated nodes
    pfix = unique(pfix,'rows') ; nfix = size(pfix,1) ;
    p = [pfix ; p];                                         % Prepend fix points

% Build the mesh   
    elmtType = pkg.geometry.mesh.elements.base.Triangle ;
    mesh = pkg.geometry.mesh.Mesh ;
    mesh.Nodes = p ;
    
% If only fixed points remains, break
    infos = {} ; % to print informations at the end
    if mesh.nNodes<=nfix ; infos{end+1} = 'out criterion: only fixed points'; end
    
% Init mesh. vizu
    meshTag = 'DistMeshPreview' ;
    delete(findobj(gca,'tag',meshTag)) ;
    if showMesh
        meshPlot = mesh.plot() ; 
        meshPlot.Tag = meshTag ; 
        meshPlot.VisibleNodes = 'all' ;
        meshPlot.BoundaryEdges.EdgeColor = 'r' ;
        meshPlot.BoundaryEdges.LineWidth = 2 ;
        meshPlot.Selected.Nodes = 1:nfix ;
        if debug ; meshPlot.UpdateOnMeshChange = true ; end
    end
    
% First mesh build
    buildMesh() ;
    
    
optimize = true ; % numel(p)~=numel(pfix) ; % Do not optimize if ther is only fixed points

if optimize

    % Initialize data
        lastPlotTime = tic ;
        lastDisplayTime = lastPlotTime ;
        count = 0 ;
        nReTri = 0 ;
    
    % OPTIMIZATION LOOP
        while optimize
            
            % Kinematic Constraints
                T = kinematicContraints ;
                nDOF = size(T,2) ;
                
            % Cost Functions
                J = sparse(nDOF,1) ;   
                H = sparse(nDOF,nDOF) ; 
                
                [rc,drc_dx] = edgeLengthCost ;
                drc_dx = drc_dx*T ;
                J = J + drc_dx'*rc ;
                H = H + drc_dx'*drc_dx ;
                
%                 [rc,drc_dx] = laplacianSmoothingCost ;
%                 drc_dx = drc_dx*T ;
%                 J = J + drc_dx'*rc ;
%                 H = H + drc_dx'*drc_dx ;
                
            % Displacement update
                if 0
                    tol = 1e-1;
                    maxIt = 10;
                    options = [] ;
                    options.type = 'nofill';
                    options.michol = 'on';
                    Lh = ichol(H);
                    [dp,~] = pcg(H,J,tol,maxIt,Lh,Lh') ;
                else
                    dp = H\J ;
                end
                dp = deltat*T*dp ;
                dp = reshape(dp,mesh.nNodes,mesh.nCoord) ;
                mesh.Nodes = mesh.Nodes - dp ;
                
            % Rebuild the mesh IF NEEDED
                dp2 = sqrt(max(sum(dp.^2,2)./fh(mesh.Nodes).^2)) ;
                if dp2>ttol
                    buildMesh() ;
                    nReTri = nReTri+1 ;
                end

            % 8. Termination criteria:
                count = count+1 ; 
                if dp2<dptol ; infos{end+1} = 'out criterion: |dP|<tol' ; break; end
                %if ~isvalid(meshPlot) ; infos{end+1} = 'out criterion: TriMesh not valid' ; break ; end
                if count>=maxCount ; infos{end+1} = 'out criterion: iteration count' ; break ; end 
                %if this.Stop.Value ; infos{end+1} = 'out criterion: user stop' ; break ; end 

            % Console display
                if displayFreq && toc(lastDisplayTime)>1/displayFreq
                    disp(['DistMesh:' ...
                            , ' it: ' , num2str(count),'/',num2str(maxCount) ...
                            , ' | Nodes: ' , num2str(mesh.nNodes) ...
                            , ' | Triangles: ' , num2str(mesh.nElems) ...
                            , ' | Re-Tri: ' , num2str(nReTri) ...
                            , ' | dP: ' , num2str(dp2,3),'/',num2str(dptol,3) ...
                    ]) ;
                end
                
            % 5. Graphical output of the current mesh
                if showMesh && toc(lastPlotTime)>1/plotFreq
                    meshPlot.Mesh = mesh ;
                    drawnow ;
                    lastPlotTime = tic ;
                end


        end
    % END OF THE OPTIMIZATION LOOP
    
    % Gather infos
        if displayFreq
            infos{end+1} = [num2str(count),' iterations'] ;
            infos{end+1} = ['last dP: ',num2str(dp2)] ;
            infos{end+1} = [num2str(nReTri),' re-triangulations'] ;
        end
    % Final Iteration
        if showMesh ; meshPlot.update ; end
        
end

% Sort elements
    mesh.sortElems ;

% DISPLAY INFOS
    if displayFreq
        infos{end+1} = [num2str(mesh.nNodes),' nodes'] ;
        infos{end+1} = [num2str(mesh.nElems),' elements'] ;
        disp(strjoin(infos,newline)) ;
    end
    
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function buildMesh()
        % Init node modifications
            validNodes = true(mesh.nNodes,1) ;
            newNodes = [] ;
        % Cull out-of-boundaries
            validNodes = validNodes & fd(mesh.Nodes)<p_dmax.*fh(mesh.Nodes) ;
        % Add/Remove nodes regarding the relative edge length
            if mesh.nEdges>0 && (tooShortThrs>eps || tooLongThrs<1/eps)
            % Relative edge length
                xe = reshape(mesh.Nodes(mesh.Edges.NodeIdx,:),mesh.nEdges,2,mesh.nCoord) ;
                Le2 = sum(diff(xe,1,2).^2,3) ; % current edge lengths
                xe = mesh.centroid(mesh.Edges) ; % edge centroids
                relLength2 = Le2./fh(xe).^2 ; % relative edge lengths
            % Cull nodes associated to too short edges
                tooShort = relLength2<tooShortThrs^2 ;
                if any(tooShort)
                    e2n = mesh.edge2node ;
                    e2n(:,~tooShort) = 0 ; % keep only too short edges
                    e2n(1:nfix,:) = 0 ; % remove fixed nodes from the analysis
                    [nn,ee] = find(e2n) ; % find free nodes attached to short edges
                    [~,is] = sort(relLength2(ee),'ascend') ; % sort by edge shortness
                    nn = nn(is) ; ee = ee(is) ;
                    [~,un] = unique(ee,'stable') ; % delete only one node by short edge
                    validNodes(nn(un)) = false ;
                end
            % Add nodes where edges are too long
            % New nodes are the centroid of big triangles
            % Because splitting edges introduces bad quality triangles
                tooLong = relLength2>tooLongThrs^2 ;
                if any(tooLong)
                % Search for triangles with all edges too long
                    ele2edg = mesh.elem2edge ;
                    tooBig = sum(ele2edg(tooLong,:),1)'==mesh.Elems.nEdges ;
                    midPt = mesh.centroid(mesh.Elems.subpart(tooBig)) ;
                    midPt = midPt(fd(midPt)<p_dmax.*fh(midPt),:) ;
                    newNodes = [newNodes ; midPt] ;
                end
            end
        % Keep fixed nodes
            validNodes(1:nfix) = true ;
        % Set the new mesh nodes
            mesh.Nodes = [mesh.Nodes(validNodes,:) ; newNodes] ;
        % Delaunay triangulation with the new nodes
            if ~all(validNodes) || ~isempty(newNodes) || mesh.nElems==0 || any(mesh.detJacobian<=0)
                tri = delaunay(mesh.Nodes) ;
                idx = padarray(tri,[0 1],1,'pre') ;
                elems = pkg.geometry.mesh.elements.ElementTable('Types',elmtType,'Indices',idx) ;
                mesh.Elems = elems ; 
            end
        % Mesh elements to keep 
            validElems = true(mesh.nElems,1) ;
        % Cull features outside the boundary
            Xt = mesh.centroid ;
            rdXt = fd(Xt)./fh(Xt) ; % relative signed distance of the triangle centroid
            validElems = validElems & rdXt<=t_dmax ;
        % Remove triangles with a quality < qmin
        % Quality is the ratio between outside and inside circle radius*2
            if qmin
                xe = reshape(mesh.Nodes(mesh.Elems.NodeIdx,:),mesh.nElems,3,mesh.nCoord) ;
                Le = sqrt(sum((xe-xe(:,[2 3 1],:)).^2,3)) ;
                a = Le(:,1) ; b = Le(:,2) ; c = Le(:,3) ;
                q = (b+c-a).*(c+a-b).*(a+b-c)./(a.*b.*c) ;
                validElems = validElems & q>=qmin ;
            end
        % Verify that no fixed point is lost
            if ~all(validElems)
                e2n = mesh.elem2node ;
                validNodes = any(e2n(:,validElems),2) ;
                cullFix = find(~validNodes(1:nfix)) ; % fixed points that should have been removed
            % If needed, reactivate the attached triangle with the minimum signed distance
                if ~isempty(cullFix)
                    d = e2n(cullFix,:).*rdXt(:)' ; % distance of attached triangles centroids
                    [~,tmin] = min(d,[],2) ;
                    validElems(tmin) = true ;
                end
            end
        % Cull triangles
            if ~all(validElems)
                mesh.Elems = mesh.Elems.subpart(validElems) ;
            % Cull unused nodes
                validNodes = ismember(1:mesh.nNodes,mesh.Elems.NodeIdx(:)) ;
                validNodes(1:nfix) = true ; % keep fixed points
                if ~all(validNodes)
                    mesh.removeNodes(~validNodes) ;
                end
            end
        % Bring the boundary points on boundary edges
            if bndCons
                ipout = mesh.boundaryNodes ;
                ipout(1:nfix) = false ; % do not move fixed points
                if any(ipout) 
                % Local gradient
                    dout = fd(mesh.Nodes(ipout,:)) ;
                    dgrad = gradfd(mesh.Nodes(ipout,:)) ;
                % Correction
                    mesh.Nodes(ipout,:) = mesh.Nodes(ipout,:)-dout.*dgrad ;
                end
            end
    end


    function [r,dr_dx] = edgeLengthCost
    % Edge Length difference
    % phi(x) = || Le(x) - Fs*fh(xm) ||²
    % Le(x) is the current edge length
    % xm is the edge centroid
        xe = mesh.Edges.dataAtIndices(mesh.Nodes) ;
        % Current Lengths
            dx = reshape(diff(xe,1,2),mesh.nEdges,mesh.nCoord)+eps ; 
            L = sqrt(sum(dx.^2,2)) ;
        % Target lengths (at middle of edges)
            xmm = reshape(mean(xe,2),mesh.nEdges,mesh.nCoord) ;
            Lt = Fscale*fh(xmm) ;
        % Residual
            r = L(:)-Lt(:) ;
        % Gradient
            ii = repmat((1:mesh.nEdges)',[1 4]) ; % [nEdges 4]
            jj = repmat(double(mesh.Edges.NodeIdx),[1 2]) + [0 0 1 1]*mesh.nNodes ; % [nEdges 4]
            dL_dx = reshape(dx,mesh.nEdges,1,mesh.nCoord).*([-1 1]./L(:)) ; % [nEdges 2 2]
            dr_dx = sparse(ii(:),jj(:),dL_dx(:),mesh.nEdges,2*mesh.nNodes) ;
    end


    function [r,dr_dx] = laplacianSmoothingCost
    % Laplacian smoothing
    % phi(x) = || x - Xc(x) ||²
    % Xc(x) is the mean of the attached triangle's barycenters
    % weighted by the element size (area)
        e2n = mesh.elem2node ;
        Mn = (e2n./sum(e2n,1))' ; % mean over nodes
        W = e2n*sparse(1:mesh.nElems,1:mesh.nElems,mesh.elemSize) ; % weights
        Me = W./sum(W,2) ; % weighted mean over elements
        Xc = Me*Mn*mesh.Nodes ;
        r = mesh.Nodes(:)-Xc(:) ;
        dr_dx = speye(2*mesh.nNodes) ;
    end

    function T = kinematicContraints
        ii = [] ; jj = [] ; vv = [] ; 
        % Get boundary nodes to constrain
            ibnd = [] ;
            if bndCons
                ibnd = find(mesh.boundaryNodes & (1:mesh.nNodes)'>nfix) ;
            end
        % Interior points only are free to move
            ifree = setdiff(1:mesh.nNodes,[1:nfix,ibnd(:)']) ;
            ii = [ii ifree ifree+mesh.nNodes] ;
            jj = [jj numel(jj)+(1:2*numel(ifree))] ;
            vv = [vv ones(1,2*numel(ifree))] ;
        % Boundary points: can only slide along the boundary tangent
            if ~isempty(ibnd)
                dgrad = gradfd(mesh.Nodes(ibnd,:)) ;
                ii = [ii ibnd(:)' ibnd(:)'+mesh.nNodes] ;
                jj = [jj repmat(numel(jj)+(1:numel(ibnd)),[1 2])] ;
                tang = flip(dgrad,2).*[-1 1] ; % [-dgrady dgrax]
                vv = [vv reshape(tang,1,[])] ;
            end
        % Transfer matrix
            T = sparse(ii(:),jj(:),vv(:),2*mesh.nNodes,max(jj(:))) ;
    end

    function dgrad = gradfd(p)
    % Return the distance function gradient at given points p
        dgrad = lvlst.gradient(p,true,geps) ;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF THE CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

