function d = toBox(P,CENTER,SIDES)
%TOBOX SIGNED Distance from a list of points P to a box
% inputs:
%   - P [nP nCoord]
%   - CENTER [1 nCoord 1] center of the rectangle
%   - SIDES [1 nCoord 1] side dimensions of the rectangle
% output: 
%   - d [nP 1 1] signed distance (d<0 inside, d>0 outside)

    if isempty(P) ; d = [] ; return ; end

    ds = abs(P-CENTER)-SIDES/2 ; % 1D Distance to box sides
    
    % INSIDE: closest side
    d = max(ds,[],2) ; % Closest side fo
    
    % OUTSIDE: edge/vertice distance
    isout = any(ds>0,2) ;
    ds(ds<=0) = 0 ;
    d(isout) = sqrt(sum(ds(isout,:).^2,2)) ;

end

function test
%% TEST THE DISTANCE
    nCoord = 2 ; nP = 1000000 ;
    CENTER = 2*rand(1,nCoord)-1 ; 
    SIDES = rand(1,nCoord) ;
    
    P = 2*(rand(nP,nCoord)-0.5)*2 ;
 %
    profile on
    tic ; d = pkg.geometry.distance.point.toBox(P,CENTER,SIDES) ; toc
    profile off
    %d(d<0) = NaN ;
    
    cla ; axis equal ; axis tight
    P3D = padarray(P,[0 max(0,3-nCoord)],0,'post') ;
    s = scatter3(P3D(:,1),P3D(:,2),P3D(:,3)+d*0,5,d,'.') ;
    s.MarkerEdgeAlpha = .01 ;
    %s.MarkerEdgeAlpha = 'flat' ; s.AlphaData = double(abs(d)<=.025) ;
    colormap(jet(50))

end



