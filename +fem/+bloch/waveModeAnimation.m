function waveModeAnimation(mesh,Kdir,U,plothandle,extrude)
% INTERACTIVE ANIMATION OF BLOCH WAVE MODES
    tag = 'waveAnim' ;
    amp = 5*15/1000./norm(range(mesh.Nodes,1)) ;
    timerPeriod = .05 ;
    animFreq = .5  ;
    Le = 3.5*norm(range(mesh.boundingBox,1)) ;
    if nargin>=5 && ~any(extrude) ; Le = 0 ; end
    de = median(mesh.elemSize(mesh.Edges)) ;

    % Delete all previous objects
    delete(timerfindall('tag',tag)) ;
    delete(findall(0,'tag',tag)) ;
    
    % Create the cursor
    fig = gcf ; ax = gca ;
    curs = patch('vertices',NaN(1,3)...
                ,'faces',[1 1 1]...
                ,'marker','.'...
                ,'EdgeColor','r'...
                ,'MarkerSize',20 ...
                ,'tag',tag) ;
    
    % Data points
    pts = [plothandle.XData(:) plothandle.YData(:) plothandle.ZData(:)] ;
                            
    % Animation figure
    figure('tag',tag) ;
    axis equal tight off ; view([30 30]) ;
    
    % Display a 3D mesh
    mesh3D = copy(mesh) ;
    if Le>0
        for cc = mesh.nCoord+1:3
            mesh3D = mesh3D.extrude(Le*full(sparse(1,mesh3D.nCoord+1,1,1,3)),ceil(Le/de)) ;
        end
    end
    mesh3D.nCoord = 3 ;
    pl0 = plot(mesh3D,'VisibleFaces','none','EdgeWidth',.05,'VisibleEdges','none') ;
    pl = plot(mesh3D,'VisibleEdges','none') ;
    if Le==0
        pl.VisibleEdges = 'all' ;
    end
    
    % 3D mode
    mode3D = @(pp)repmat(amp*U(:,:,pp)/max(reshape(abs(U(:,:,pp)),[],1)),[mesh3D.nNodes/mesh.nNodes 1]).*exp(-1i*sum(Kdir(pp,:).*mesh3D.Nodes,2,'omitnan')) ;
    curs.UserData = mode3D(1) ;
    
    % Cursor function
    fig.WindowButtonMotionFcn = @(src,evt)cellfun(@(c)c(src,evt),{...
                    @(src,evt)set(curs,'vertices',pts(pkg.ui.closestToMouse(pts,ax),:)) ...
                    @(src,evt)set(curs,'UserData',mode3D(pkg.ui.closestToMouse(pts,ax))) ...
                                },'uni',false);

    % Timer that makes the wave phase turn                        
    startTime = tic ;
    defShape = @()real(curs.UserData*exp(2i*pi*toc(startTime)*animFreq)) ;
    timerfcn = @(src,evt)cellfun(@(c)c(src,evt),{...
                    @(src,evt)set(pl,'Deformation',defShape(),'CData',sqrt(sum(defShape().^2,2))) ...
                },'uni',false);
    ti = timer('Period',timerPeriod,'ExecutionMode','FixedRate','TimerFcn',timerfcn,'tag',tag) ;
    start(ti)

    % Delete timer on figure(s) close
    addlistener(pl.Parent,'ObjectBeingDestroyed',@(src,evt)delete(ti)) ;
    addlistener(pl.Parent,'ObjectBeingDestroyed',@(src,evt)set(fig,'WindowButtonMotionFcn',[])) ;
%     addlistener(curs.Parent,'ObjectBeingDestroyed',@(src,evt)delete(ti)) ;
%     addlistener(curs.Parent,'ObjectBeingDestroyed',@(src,evt)set(fig,'WindowButtonMotionFcn',[])) ;

end