function mesh = distMesh(lvlst,varargin)
%DISTMESH build a distance mesh from a pkg.geometry.levelset.LevelSet
% mesh = distMesh(lvlst) automatic meshing
% mesh = distMesh(lvlst,h0) uniform edge length h0 (scalar)
% mesh = distMesh(lvlst,'Name','Value',...) custom parameters

% Process input
    if ~isa(lvlst,'pkg.geometry.levelset.LevelSet')
        error('The first argument MUST be a pkg.geometry.levelset.LevelSet object') ;
    end
    if mod(nargin,2)==0 % given edge length
        varargin = [{'h0'} varargin] ;
    end
    args = struct(varargin{:}) ; 

% Default input arguments 
    % Distance function
        if isfield(args,'fd') ; fd = args.fd ;
        else ; fd = @lvlst.Function ;
        end
    % Mesh boundary only ?
        if isfield(args,'bnd_only') ; bnd_only = args.bnd_only ;
        else ; bnd_only = false ;
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
            else ; p0 = lvlst.populate(h0,'iso',fh,bnd_only) ;
            end
        % Fixed nodes
            if isfield(args,'pfix') ; pfix = args.pfix ;
            else ; pfix =  lvlst.discretizeEdges(fh) ... % variable-spaced contour points
                        ... lvlst.discretizeEdges(h0) ... % uniform spaced contour points
                        ... lvlst.Kinks ... % only kink points
                        ... [] ... % no points
                        ;
            end
    % RELAXATION PARAMETERS
        % Max number of iterations
            if isfield(args,'maxCount') ; maxCount = args.maxCount ;
            else ; maxCount = 100 ;
            end
        % Maximum allowed displacement (relative to local edge length
            if isfield(args,'maxdp') ; maxdp = args.maxdp;
            else ; maxdp = 1.0 ;
            end
        % Node displacement tolerance (convergence criterion, relative to local edge length)
            if isfield(args,'dptol') ; dptol = args.dptol;
            else ; dptol = 0.01 ;
            end
        % Re-triangulation tolerance (relative to local edge length)
            if isfield(args,'ttol') ; ttol = args.ttol ;
            else ; ttol = 0.2 ;
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
        % Global displacement regularization (try to avoid solid motion of symmetric shapes)
            if isfield(args,'reg_global') ; reg_global = args.reg_global ;
            else ; reg_global = 1e-3 ;
            end
        % Maximum distance function allowed for Nodes (relative to local edge length)
            if isfield(args,'p_dmax') ; p_dmax = args.p_dmax ;
            else ; p_dmax =  0.5 ; 0.001 ;
            end
        % Maximum distance function allowed for Element centroids (relative to local edge length)
            if isfield(args,'t_dmax') ; t_dmax = args.t_dmax ;
            else 
                if bnd_only ; t_dmax = -0.02 ;
                else ; t_dmax = -0.2 ;
                end
            end
        % Too SHORT edges threshold (rel. to local edge length)
            if isfield(args,'tooShortThrs') ; tooShortThrs = args.tooShortThrs ;
            else ; tooShortThrs = 0.75 ;
            end
        % Too LONG edges threshold (rel. to local edge length)
            if isfield(args,'tooLongThrs') ; tooLongThrs = args.tooLongThrs ;
            else ; tooLongThrs = 1.25 ;
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
    
% Delete too close nodes (e.g if p0 has close duplicates with pfix)
    [~,ia] = uniquetol(p,lvlst.defaultTolerance,'ByRows',true,'DataScale',1,'OutputAllIndices',true) ;
    ia = cellfun(@min,ia) ; % take the first duplicatd point, avoiding the deletion of any pfix
    p = p(sort(ia),:) ;

% Build the mesh   
    switch lvlst.nCoord
        case 3
            elmtType = pkg.geometry.mesh.elements.base.Tetrahedron ;
        otherwise
            elmtType = pkg.geometry.mesh.elements.base.Triangle ;
    end
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
        lastBuildMeshNodes = mesh.Nodes ;
    
    % OPTIMIZATION LOOP
        while optimize
            
            % Kinematic Constraints
                T = kinematicContraints ;
                nDOF = size(T,2) ;
                if 1 % visualize degrees of freedom
%                     dofTag = [meshTag '_DOF'] ;
%                     delete(findobj(gca,'tag',dofTag)) ;
%                     [nnc,ddd,vvv] = find(T) ; % [dof node*coord value]
%                     [nnn,ccc] = ind2sub(size(mesh.Nodes),nnc) ; % [node coord]
%                     ndof = accumarray(ddd,nnn,[nDOF 1],@min) ; % DOF corresponding node
%                     vdof = full(sparse(ddd,ccc,vvv)) ; % DOF vector
%                     xdof = mesh.Nodes(ndof,:) ;
%                     xdof(:,end+1:3) = 0 ; % force 3D coords
%                     vdof(:,end+1:3) = 0 ; % force 3D coords
%                     quiver3(xdof(:,1),xdof(:,2),xdof(:,3) ...
%                             ,vdof(:,1),vdof(:,2),vdof(:,3) ...
%                             ,'tag',dofTag)
                end
                
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

            % Global displacement regularization (try to avoid solid motion)
                if reg_global
                    Hd = T'*T ;
                    H = H + reg_global*Hd*sum(abs(H(:)))./sum(abs(Hd(:))) ;
                end
                
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
                
            % Maximum displacement reached ?
                if ~isinf(maxdp)
                    dpmax = sqrt(max(sum(dp.^2,2)./fh(mesh.Nodes).^2)) ;
                    if dpmax>maxdp ; dp = maxdp*(dp./dpmax) ; end
                end
                
            % If needed, correct points on boundary
                if bnd_only
                    ptemp = mesh.Nodes - dp ;
                    ptemp = lvlst.toBoundary(ptemp) ;
                    dp = mesh.Nodes - ptemp ;
                end
                
            % Apply to the mesh
                mesh.Nodes = mesh.Nodes - dp ;
                dp2 = max(sum(dp.^2,2)./fh(mesh.Nodes).^2) ;
                
            % Rebuild the mesh IF NEEDED
                dpMesh2 = max(sum((mesh.Nodes-lastBuildMeshNodes).^2,2)./fh(mesh.Nodes).^2) ;
                if dpMesh2>ttol^2
                    buildMesh() ;
                    lastBuildMeshNodes = mesh.Nodes ;
                    nReTri = nReTri+1 ;
                end

            % 8. Termination criteria:
                count = count+1 ; 
                if dp2<dptol^2 ; infos{end+1} = 'out criterion: |dP|<tol' ; break; end
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
                            , ' | dP: ' , num2str(sqrt(dp2),3),'/',num2str(dptol,3) ...
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
            infos{end+1} = ['last dP: ',num2str(sqrt(dp2))] ;
            infos{end+1} = [num2str(nReTri),' re-triangulations'] ;
        end
    % Final Iteration
        if showMesh ; meshPlot.update ; end
        
end

% Sort elements
    if mesh.nCoord==2 ; mesh.sortElems ; end

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
            if bnd_only ; newNodes = lvlst.toBoundary(newNodes) ; end
            mesh.Nodes = [mesh.Nodes(validNodes,:) ; newNodes] ;
        % Delaunay (tri.tet)angulation with the new nodes
            remesh = ~all(validNodes) || ~isempty(newNodes) || mesh.nElems==0 || any(mesh.detJacobian<=0) ;
            if remesh
%                 tri = delaunayn(mesh.Nodes,{'QJ','Qt','Qbb','Qc'}) ;
                tri = delaunay(mesh.Nodes) ;
                idx = padarray(tri,[0 1],1,'pre') ;
                elems = pkg.geometry.mesh.elements.ElementTable('Types',elmtType,'Indices',idx) ;
                mesh.Elems = elems ; 
            end
        % Mesh elements to keep 
            validElems = true(mesh.nElems,1) ;
        % Cull features outside the boundary
            if ~bnd_only || remesh
                Xt = mesh.centroid ;
                rdXt = fd(Xt)./fh(Xt) ; % relative signed distance of the triangle centroid
                validElems = validElems & rdXt<=t_dmax ;
            end
        % If only the boundary has to be meshed
            if bnd_only
                mesh.Elems = mesh.Elems.subpart(validElems) ; 
                switch mesh.nCoord
                    case 2
                        mesh.Elems = mesh.Edges.subpart(mesh.boundaryEdges) ;
                    case 3
                        mesh.Elems = mesh.Faces.subpart(mesh.outerFaces) ;
                end
            % Reset valid elements & centroids 
                validElems = true(mesh.nElems,1) ;
                Xt = mesh.centroid ;
                rdXt = fd(Xt)./fh(Xt) ; % relative signed distance of the triangle centroid
            end
        % Remove simplices with a quality < qmin
        % Quality is the ratio between inside and outside circle radius
            if qmin
                q = mesh.elemQuality ;
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
        % Cull Elements
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
                switch mesh.nCoord
                    case 3 ; ipout = mesh.outerNodes ;
                    otherwise ; ipout = mesh.boundaryNodes ;
                end
                ipout(1:nfix) = false ; % do not move fixed points
                if any(ipout) 
                % Correction
                    mesh.Nodes(ipout,:) = lvlst.toBoundary(mesh.Nodes(ipout,:)) ;
                end
            end
    end


    function [r,dr_dx] = edgeLengthCost
    % Edge Length difference
    % phi(x) = || Le(x) - Fs*fh(xm) ||�
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
            ii = repmat((1:mesh.nEdges)',[1 2*mesh.nCoord]) ; % [nEdges 2*mesh.nCoord]
            jj = repmat(double(mesh.Edges.NodeIdx),[1 mesh.nCoord]) + repelem(0:mesh.nCoord-1,2)*mesh.nNodes ; % [nEdges 2*mesh.nCoord]
            dL_dx = reshape(dx,mesh.nEdges,1,mesh.nCoord).*([-1 1]./L(:)) ; % [nEdges 2 nCoord]
            dr_dx = sparse(ii(:),jj(:),dL_dx(:),mesh.nEdges,mesh.nCoord*mesh.nNodes) ;
    end


    function [r,dr_dx] = laplacianSmoothingCost
    % Laplacian smoothing
    % phi(x) = || x - Xc(x) ||�
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
                switch mesh.nCoord
                    case 3 ; ibnd = find(mesh.outerNodes & (1:mesh.nNodes)'>nfix) ;
                    case 2 ; ibnd = find(mesh.boundaryNodes & (1:mesh.nNodes)'>nfix) ;
                end
            end
            ibnd = ibnd(:)' ; % row vector
        % Interior points only are free to move
            ifree = setdiff(1:mesh.nNodes,[1:nfix,ibnd]) ;
            ii = [ii ifree ifree+mesh.nNodes] ;
            jj = [jj numel(jj)+(1:2*numel(ifree))] ;
            vv = [vv ones(1,2*numel(ifree))] ;
        % Boundary points: can only slide along the boundary tangent (curve/plane)
            if ~isempty(ibnd)
                pb = mesh.Nodes(ibnd,:) ;
                dgrad = gradfd(pb) ; % lvlst gradient == normal to the contour
                switch mesh.nCoord
                    case 2
                        tang = flip(dgrad,2).*[-1 1] ; % [-dgrady dgrax]
                        iit = [ibnd ibnd+mesh.nNodes] ;
                        jjt = repmat(1:numel(ibnd),[1 2]) ;
                    case 3
                        tang = cellfun(@null,num2cell(dgrad,2),'uni',false) ; % tangent space: nP cell array of [3 2]==[nCoord nTangents]
                        tang = cat(3,tang{:}) ; % [3==nCoord 2=nTangents nP]
                        iit = mesh.nNodes*(0:2)' + [0 0] + reshape(ibnd,1,1,[]) ;  % [3==nCoord 2=nTangents nP]
                        jjt = [0;0;0] + numel(ibnd)*[0 1] + reshape(1:numel(ibnd),1,1,[]) ;  % [3==nCoord 2=nTangents nP]
                end
                ii = [ii iit(:)'] ;
                jj = [jj numel(jj)+jjt(:)'] ;
                vv = [vv tang(:)'] ;
            end
        % Transfer matrix
            T = sparse(ii(:),jj(:),vv(:),mesh.nCoord*mesh.nNodes,max(jj(:))) ;
    end

    function dgrad = gradfd(p)
    % Return the distance function gradient at given points p
        dgrad = lvlst.gradient(p,true,geps) ;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF THE CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

