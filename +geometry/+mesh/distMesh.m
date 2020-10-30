function mesh = distMesh(lvlst,varargin)
%DISTMESH build a distance mesh from a pkg.geometry.levelset.LevelSet

if ~isa(lvlst,'pkg.geometry.levelset.LevelSet')
    error('The first argument MUST be a pkg.geometry.levelset.LevelSet object') ;
end

% Default input arguments
    fd = @lvlst.Function ; % signed distance function
    h0 = norm(range(lvlst.BoundingBox,1))/50 ; % initial edge length 
    fh = @(p)h0*(ones(size(p,1),1)) ; %+3*abs(p(:,2))) ; % space-dependent edge length
    pfix = ... % fixed points
             lvlst.discretizeContour(fh) ... % variable-spaced contour points
            ... lvlst.discretizeContour(h0) ... % uniform spaced contour points
            ... lvlst.Kinks ... % only kink points
            ... [] ... % no points
        ;
    initialDistribution = 'iso' ;
    ax = gca ; % axes where to plot the mesh update
    showMesh = ~isempty(ax) ;
    debug = false && showMesh ; %link the meshplot to any mesh change

% Parameters for convergence of the mesh
    maxCount = 100 ; % max number of iterations
    dptol = 0.01 ; % point displacement tolerance (convergence criterion, relative to mean density)
    ttol = 0.01 ; % re-triangulation tolerance
    qmin = 0*0.3 ; % Minimum triangle quality
    Fscale = 1.2 ; % spring L0 relative to bar length
    deltat = .3 ; % explicit time increment
    p_dmax = 0.5*1e-0*h0; % maximum distance function allowed (remove pts)
    t_dmax = -0.2*h0 ; % maximum distance function allowed (remove triangles)
    geps = 1e-2*h0 ; %*sqrt(eps) % discrete gradient estim. step
    tooShortThrs = 0.75 ;
    tooLongThrs = 1.25 ;
    plotFreq = Inf ; % update plot frequency
    meshTag = 'DistMeshPreview' ;

% Create initial distribution in bounding box (equilateral triangles)
    p = lvlst.populate(h0,initialDistribution,fh) ;

% Remove points outside the geometry
    p = p(fd(p)<p_dmax,:) ;

% Apply the rejection method
%     r0=1./fh(p).^2 ;                % Probability to keep point
%     pkeep = rand(size(p,1),1)<r0./max(r0) ;
%     p = p(pkeep,:) ;                % Rejection method
    
% Append fixed nodes
    if ~isempty(pfix), p=setdiff(p,pfix,'rows'); end     % Remove duplicated nodes
    pfix = unique(pfix,'rows'); nfix=size(pfix,1) ;
    p = [pfix ; p];                                         % Prepend fix points

% Build the mesh   
    elmtType = pkg.geometry.mesh.elements.base.Triangle ;
    mesh = pkg.geometry.mesh.Mesh ;
    
% If only fixed points remains, break
    infos = {} ; % to print informations at the end
    if mesh.nNodes<=nfix ; infos{end+1} = 'out criterion: only fixed points'; end
    
% Init mesh. vizu
    delete(findobj(ax,'tag',meshTag)) ;
    if showMesh
        meshPlot = mesh.plot() ; 
        meshPlot.Tag = meshTag ; 
        meshPlot.VisibleNodes = 'all' ;
        meshPlot.BoundaryEdges.EdgeColor = 'r' ;
        meshPlot.BoundaryEdges.LineWidth = 2 ;
        if debug ; meshPlot.UpdateOnMeshChange = true ; end
    end
    
% First mesh build
    buildMesh(p) ;
    
    
optimize = true ; % numel(p)~=numel(pfix) ; % Do not optimize if ther is only fixed points

if optimize

    % Initialize data
        lastPlotTime = tic ;
        count = 0 ;
        nReTri = 0 ;
    
    % OPTIMIZATION LOOP
        while optimize
                
            % Cost Functions
                J = sparse(2*mesh.nNodes,1) ;   
                H = sparse(2*mesh.nNodes,2*mesh.nNodes) ; 
                
                [rc,Jc,Hc] = edgeLengthCost ;
                J = J + Jc ;
                H = H + Hc ;
                
%                 [rc,Jc,Hc] = laplacianSmoothingCost ;
%                 J = J + 0.1*Jc ;
%                 H = H + 0.1*Hc ;
            
            % Kinematic Constraints
                T = kinematicContraints ;
                J = T'*J ;
                H = T'*H*T ;
                
            % Displacement update
                dp = deltat*T*(H\J) ;
                dp = reshape(dp,mesh.nNodes,mesh.nCoord) ;
                mesh.Nodes = mesh.Nodes - dp ;
                
            % Rebuild the mesh IF NEEDED
                dp2 = sqrt(max(sum(dp.^2,2)./fh(mesh.Nodes).^2)) ;
                if dp2>ttol
                    buildMesh(mesh.Nodes) ;
                    nReTri = nReTri+1 ;
                end

            % 8. Termination criteria:
                count = count+1 ; 
                if dp2<dptol ; infos{end+1} = 'out criterion: |dP|<tol' ; break; end
                %if ~isvalid(meshPlot) ; infos{end+1} = 'out criterion: TriMesh not valid' ; break ; end
                if count>=maxCount ; infos{end+1} = 'out criterion: iteration count' ; break ; end 
                %if this.Stop.Value ; infos{end+1} = 'out criterion: user stop' ; break ; end 

            % 
            
            
            % 5. Graphical output of the current mesh
                if toc(lastPlotTime)>1/plotFreq
                    if showMesh ; meshPlot.Mesh = mesh ; end
                    disp(['DistMesh:' ...
                            , ' it: ' , num2str(count),'/',num2str(maxCount) ...
                            , ' | Nodes: ' , num2str(mesh.nNodes) ...
                            , ' | Triangles: ' , num2str(mesh.nElems) ...
                            , ' | Re-Tri: ' , num2str(nReTri) ...
                            , ' | dP: ' , num2str(dp2,3),'/',num2str(dptol,3) ...
                    ]) ;
                    drawnow ;
                    lastPlotTime = tic ;
                end


        end
    % END OF THE OPTIMIZATION LOOP
    
    % Gather infos
        infos{end+1} = [num2str(count),' iterations'] ;
        infos{end+1} = ['last dP: ',num2str(dp2)] ;
        %infos{end+1} = [num2str(nReTri),' re-triangulations'] ;
    % Clean up and plot final mesh
        %[p,t,pix] = fixmesh(p,t) ;
    % Final Iteration
        if showMesh ; meshPlot.update ; end
        
end

% Sort elements
    mesh.sortElems ;

% Format the output
    infos{end+1} = [num2str(mesh.nNodes),' nodes'] ;
    infos{end+1} = [num2str(mesh.nElems),' elements'] ;
    
% DISPLAY INFOS
    disp(strjoin(infos,newline)) ;
    
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function buildMesh(p)
        % First delaunay triangulation for node deletion
            mesh.Nodes = p ;
            idx = padarray(delaunay(p),[0 1],1,'pre') ;
            elems = pkg.geometry.mesh.elements.ElementTable('Types',elmtType,'Indices',idx) ;
            mesh.Elems = elems ; 
        % Valid nodes
            validNodes = true(mesh.nNodes,1) ;
            newNodes = [] ;
        % Cull out-of-boundaries
            validNodes = validNodes & fd(mesh.Nodes)<p_dmax ;
        % Cull nodes associated to too short edges
            Le = mesh.elemSize(mesh.Edges) ; % current edge lengths
            xe = mesh.centroid(mesh.Edges) ; % edge centroids
            relLength = Le./fh(xe) ; % relative edge lengths
            tooShort = relLength<tooShortThrs ;
            if any(tooShort)
                e2n = mesh.edge2node ;
                e2n = e2n.*tooShort(:)' ; % keep only too short edges
                e2n = e2n.*~sparse(1:nfix,1,true,mesh.nNodes,1) ; % remove fixed nodes from the analysis
                [nn,ee] = find(e2n) ; % find free nodes attached to short edges
                [~,is] = sort(relLength(ee),'ascend') ; % sort by edge shortness
                nn = nn(is) ; ee = ee(is) ;
                [~,un] = unique(ee,'stable') ; % delete only one node by short edge
                validNodes(nn(un)) = false ;
            end
        % Keep fixed nodes
            validNodes(1:nfix) = true ;
        % Add nodes where edges are too long
        % New nodes are the centroid of big triangles
        % Because splitting edges introduces bad quality triangles
            tooLong = relLength>tooLongThrs ;
            if any(tooLong)
            % Search for triangles with all edges too long
                ele2edg = mesh.elem2edge ;
                tooBig = sum(ele2edg(tooLong,:),1)'==mesh.Elems.nEdges ;
                midPt = mesh.centroid(mesh.Elems.subpart(tooBig)) ;
                midPt = midPt(fd(midPt)<p_dmax,:) ;
                newNodes = [newNodes ; midPt] ;
            end
        % Set new Nodes
            mesh.Nodes = [mesh.Nodes(validNodes,:) ; newNodes];
        % Re-triangulate
            elems.Indices = padarray(delaunay(mesh.Nodes),[0 1],1,'pre') ;
            mesh.Elems = elems ;
        % Mesh features to keep 
            validElems = true(mesh.nElems,1) ;
        % Cull features outside the boundary
            validElems = validElems & fd(mesh.centroid)<=t_dmax ;
        % Remove triangles with a quality < qmin
        % Quality is the ratio between outside and inside circle radius*2
            if qmin
                xe = mesh.Elems.dataAtIndices(mesh.Nodes) ; % [nElems 3 nCoord]
                Le = sqrt(sum((xe-xe(:,[2 3 1],:)).^2,3)) ;
                a = Le(:,1) ; b = Le(:,2) ; c = Le(:,3) ;
                q = (b+c-a).*(c+a-b).*(a+b-c)./(a.*b.*c) ;
                validElems = validElems & q>=qmin ;
            end
        % Keep triangles attached to fixed points
            %valid = valid | logical(mesh.elem2node'*sparse(1:nfix,ones(1,nfix),true,mesh.nNodes,1)) ;
        % Cull triangles
            mesh.Elems = mesh.Elems.subpart(validElems) ;
        % Cull unused nodes
            mesh.cullUnused ;
        % Bring the boundary points on boundary edges
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


    function [r,J,H] = edgeLengthCost
    % Edge Length difference
    % phi(x) = || Le(x) - Fs*fh(xm) ||²
    % Le(x) is the current edge length
    % xm is the edge centroid
        xe = mesh.Edges.dataAtIndices(mesh.Nodes) ;
        % Current Lengths
            dx = reshape(diff(xe,1,2),mesh.nEdges,mesh.nCoord) ; 
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
        % G-N Matrices
            J = dr_dx'*r ;
            H = dr_dx'*dr_dx ;
    end


    function [r,J,H] = laplacianSmoothingCost
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
        J = dr_dx'*r ;
        H = dr_dx'*dr_dx ;
    end

    function T = kinematicContraints
        ibnd = setdiff(find(mesh.boundaryNodes),1:nfix) ;
        ii = [] ; jj = [] ; vv = [] ; 
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
    % mean gradient over a quad
        if isempty(p) ; dgrad = [] ; return ; end
        sp = [-1 -1 ; 1 -1 ; 1 1 ; -1 1]/2 ; % [4 nCoord] shifts
        pp = permute(p,[1 3 2]) + permute(sp,[3 1 2])*geps ; % [nP 4 nCoord] all points
        d = reshape(fd(reshape(pp,[],size(p,2))),size(p,1),[]) ; % [nP 4] distance values
        dgrad = (d*sp)*(1/geps) ; % [nP 2] gradient
        normGrad = sqrt(sum(dgrad.^2,2)) ;
        dgrad = dgrad./normGrad ;
        if any(normGrad<1e-1) ; error('Vanishing gradient found !') ; end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF THE CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

