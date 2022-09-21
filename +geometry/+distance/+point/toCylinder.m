function d = toCylinder(P,RADIUS,PTS,ISINF)
%TOCYLINDER SIGNED Distance from a list of points P to a (list of) (IN)FINITE Cylinder(s)
% inputs:
%   - P [nP nCoord] list of points
%   - RADIUS [1 nCoord nC] radius of the cylinder(s)
%   - PTS [1 nCoord nC] center points of the cylinder face(s)
%   - ISINF [1 1 nC] is the cylinder infinite ?
% output: 
%   - d [nP 1 nC] signed distance (d<0 inside, d>0 outside)

    if isempty(P) ; d = [] ; return ; end
    if nargin<4 ; ISINF = false(1,1,size(RADIUS,3)) ; end
    
    % Distance to infinite cylinders
    d = pkg.geometry.distance.point.toLine(P,PTS(1,:,:),PTS(2,:,:))-RADIUS ; % [nP 1 nC]
    
    if all(ISINF) % All cylinders are infinite ?
        return ;
    end
    
    
    % Memory access optimization
    ISINF = ISINF & true(size(RADIUS)) ;
    dc = d(:,:,~ISINF) ;
    PTS = PTS(:,:,~ISINF) ; % [2 nCoord nCf]
    P1 = PTS(1,:,:) ; % [1 nCoord nCf]
    P2 = PTS(2,:,:) ; % [1 nCoord nCf]
    
    % Distance to top/bottom faces planes
    AX = P2-P1 ; % oriented cylinder axis [1 nCoord nC]
    df1 = pkg.geometry.distance.point.toPlane(P,P1,-AX) ; % [nP 1 nCf]
    df2 = pkg.geometry.distance.point.toPlane(P,P2,AX) ; % [nP 1 nCf]
    
    % Which point is inside which feature ? [infiniteCylinder, facePlane1, facePlane2]
    isin = [dc df1 df2]<=0 ; % [nP 3 nCf]
    %allout = all(~isinf,2) ; % this should be all zero !
    
    % Distance to finite cylinder: 
    % if inside infinite cylinder & all faces, maximum signed distance with these features
    allin = all(isin,2) ; % [nP 1 nCf]
    dc(allin) = max([dc(allin) df1(allin) df2(allin)],[],2) ;
    % if inside the infinite cylinder but outside a face, take the face distance
    inC_outF1 = all(isin==[1 0 1],2) ;
    dc(inC_outF1) = df1(inC_outF1) ;
    inC_outF2 = all(isin==[1 1 0],2) ;
    dc(inC_outF2) = df2(inC_outF2) ;
    % if outside the cylinder and outside ONE face, take the distance to the face's edge
    outC_outF1 = all(isin==[0 0 1],2) ;
    dc(outC_outF1) = sqrt(dc(outC_outF1).^2 + df1(outC_outF1).^2) ;
    outC_outF2 = all(isin==[0 1 0],2) ;
    dc(outC_outF2) = sqrt(dc(outC_outF2).^2 + df2(outC_outF2).^2) ;
    % if outside inifite cylinder and inside both faces, keep the distance to axis
    
    % Set
    d(:,:,~ISINF) = dc ;

end

function test
%% TEST THE DISTANCE
    nC = 10 ; nP = 1000000 ;
    R = rand(1,1,nC) ;
    PTS = rand(2,3,nC) ;
    
    P = 2*(rand(nP,3)-0.5)*2.*[1 1 1] + [.5 .5 .5] ;
    
    profile on
    tic ; d = pkg.geometry.distance.point.toCylinder(P,R,PTS,false) ; toc
    profile off
    
    cla ; axis equal
    s = scatter3(P(:,1),P(:,2),P(:,3),5,min(d,[],3),'.') ;
    %s.MarkerEdgeAlpha = .01 ;
    s.MarkerEdgeAlpha = 'flat' ; s.AlphaData = double(abs(min(d,[],3))<=.025) ;

end

