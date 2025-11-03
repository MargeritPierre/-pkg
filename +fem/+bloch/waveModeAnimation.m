function [pl] = waveModeAnimation(mesh,Kdir,U,plothandle,extrude_L_or_N,gif_on_click)
% INTERACTIVE ANIMATION OF BLOCH WAVE MODES
    tag = "waveAnim" ;
    amp = .05.*norm(range(mesh.Nodes,1)) ;
    timerPeriod = .05 ;
    animFreq = .5  ;
    
% Unit cell replication option
    if nargin>=5 && any(extrude_L_or_N)
        if islogical(extrude_L_or_N) % extrude==true
            extrude_L_or_N = 3.5*norm(range(mesh.boundingBox,1)) ;
%         elseif isinteger(extrude_L_or_N) % number of unit cells
%             extrude_L_or_N = 0 ; 
        end
    else
        extrude_L_or_N = 0 ;
    end

    % Delete all previous objects
    delete(timerfindall('tag',tag)) ;
    delete(findall(0,'tag',tag)) ;
    
    % Create the cursor
    ax = plothandle.Parent ; fig = ax.Parent ;
    curs = patch(ax,'vertices',NaN(1,3)...
                ,'faces',[1 1 1]...
                ,'marker','.'...
                ,'EdgeColor','r'...
                ,'MarkerSize',20 ...
                ,'handlevisibility','off' ...
                ,'tag',tag) ;
    
    % Data points
    pts = [plothandle.XData(:) plothandle.YData(:) plothandle.ZData(:)] ;
    
    % Build the 3D mesh
    wavemesh = copy(mesh) ;
    if any(extrude_L_or_N>0)
        if size(Kdir,2)==mesh.nCoord % architectured material, tile the unit cell
            extrude_N = uint8(extrude_L_or_N).*ones(1,mesh.nCoord,'uint8') ;
            repmesh = pkg.geometry.mesh.GridMesh(extrude_N-1) ; 
            x = permute(repmesh.Nodes,[3 2 1]).*range(mesh.boundingBox(),1) ;
            wavemesh = mesh.replicate(mesh.Nodes+x,false) ;
        else % waveguides, extrude along the remaingin directions
            de = median(mesh.elemSize(mesh.Edges)) ;
            for cc = mesh.nCoord+1:3
                wavemesh = wavemesh.extrude(extrude_L_or_N*full(sparse(1,wavemesh.nCoord+1,1,1,3)),ceil(extrude_L_or_N/de)) ;
            end
        end
    end
    % Extend coordinates if needed
        U(:,end+1:wavemesh.nCoord,:) = 0 ;
        Kdir(:,end+1:wavemesh.nCoord) = 0 ;
    
    % Normalize
    U = amp*U./max(abs(U),[],1:2) ;
                            
    % Animation figure
    fiig = findall(0,'tag',tag+"_fig") ;
    if isempty(fiig) ; fiig = figure() ; end % create a new figure only if needed
    clf(fiig,'reset') ;
    set(fiig,'tag',tag+"_fig") ;
    axis equal tight off ; 
    set(gca,'Clipping','off')
    if wavemesh.nCoord>2 ; view([30 30]) ; end
    caxis([-1 1]*amp) ;
        
    % Display
    pl0 = plot(wavemesh,'VisibleFaces','none','EdgeWidth',.05,'VisibleEdges','none') ;
    pl = plot(wavemesh,'VisibleEdges','none') ;
    if extrude_L_or_N==0
        pl.VisibleEdges = 'all' ;
    end
    
    axLims = wavemesh.boundingBox()+amp.*[-1.1;1.1] ;
    axWave = gca ;
    axWave.XTick = [] ; axWave.YTick = [] ; axWave.ZTick = [] ;
    axWave.XLim = axLims(:,1) ; axWave.YLim = axLims(:,2) ; 
    if wavemesh.nCoord>2 ; axWave.ZLim = axLims(:,3) ; end
    
    % 3D mode
    mode3D = @(pp)repmat(U(:,:,pp),[wavemesh.nNodes/mesh.nNodes 1]).*exp(-1i*sum(Kdir(pp,:).*wavemesh.Nodes,2,'omitnan')) ;
    normalize = @(u)u./max(abs(u(:))) ;
    curs.UserData = mode3D(1) ;
    
    % Cursor function
    fig.WindowButtonMotionFcn = @(src,evt)cellfun(@(c)c(src,evt),{...
                    @(src,evt)set(curs,'vertices',pts(pkg.ui.closestToMouse(pts,ax),:)) ...
                    @(src,evt)set(curs,'UserData',amp*normalize(mode3D(pkg.ui.closestToMouse(pts,ax)))) ...
                                },'uni',false);

    % Timer that makes the wave phase turn                        
    startTime = tic ;
    defShape = @(t)real(curs.UserData*exp(2i*pi*t*animFreq)) ;
    timerfcn = @(src,evt)cellfun(@(c)c(src,evt),{...
                    @(src,evt)set(pl,'Deformation',defShape(toc(startTime)),'CData',sqrt(sum(defShape(toc(startTime)).^2,2))) ...
                },'uni',false);
    ti = timer('Period',timerPeriod,'ExecutionMode','FixedSpacing','TimerFcn',timerfcn,'tag',tag) ;
    start(ti)

    % Delete timer on figure(s) close
    addlistener(pl.Parent,'ObjectBeingDestroyed',@(src,evt)delete(ti)) ;
    addlistener(pl.Parent,'ObjectBeingDestroyed',@(src,evt)set(fig,'WindowButtonMotionFcn',[])) ;
    addlistener(curs.Parent,'ObjectBeingDestroyed',@(src,evt)delete(ti)) ;
    addlistener(curs.Parent,'ObjectBeingDestroyed',@(src,evt)set(fig,'WindowButtonMotionFcn',[])) ;

    % record gif on click
    if nargin>5 && gif_on_click
        fig.WindowButtonDownFcn = @(src,evt)cellfun(@(c)c(src,evt),{...
                        @(src,evt)clickFcn ...
                    },'uni',false);
        addlistener(pl.Parent,'ObjectBeingDestroyed',@(src,evt)set(fig,'WindowButtonDownFcn',[])) ;
        frame = annotation(gcf,'rectangle',[0 0 1 1],'linewidth',5,'Color','none') ;
    end

    function clickFcn
        cc = copyobj(curs,curs.Parent) ;
        cc.Tag = 'mode' ;
        nGifs = numel(findall(cc.Parent,'tag','mode')) ;
        clr = cc.Parent.ColorOrder(mod(nGifs-1,end-1)+1,:) ;
        cc.EdgeColor = clr ;
        frame.Color = clr ;
        % txt = text(cc.XData,cc.YData,cc.ZData...
        %                 ,"\boldmath $" + num2str(nGifs) + "$" ...
        %                 ,'Color','w' ...
        %                 ,'BackgroundColor',cc.EdgeColor ...
        %                 ) ;

        stop(ti) ;
        IMG = {} ;
        for tt = 0:timerPeriod:1/animFreq-eps
            set(pl,'Deformation',defShape(tt),'CData',sqrt(sum(defShape(tt).^2,2)))
            IMG{end+1} = getframe(fiig) ;
        end
        filename = "mode_"+num2str(nGifs)+".gif" ;
        pkg.export.gif(IMG,filename,'DelayTime',timerPeriod)
        frame.Color = 'none' ;
        start(ti) ;
    end
end



