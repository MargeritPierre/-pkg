function d = toPolytope(P,O,N,approx)
%TOPOLYTOPE Signed distance from a list of points P to a CONVEX polytope
% The polytope is defined by a set of oriented planes (that form the faces
% of the polytope)
% inputs:
%   - P [nP nCoord] list of points
%   - O [nCP nCoord] oriented plane origin points
%   - N [nCP nCoord] oriented plane outgoing normals
%   - approx: do not include edge distance
% output: 
%   - signed distance d [nP 1 1] 

    if isempty(P) ; d = [] ; return ; end
    if nargin<4 ; approx = false ; end
    %nP = size(P,1) ;
    nCP = size(O,1) ; % number of "cut" planes
    nCoord = size(P,2) ;

    O = permute(O,[3 2 1]) ; % [1 nCoord nCP]
    N = permute(N,[3 2 1]) ; % [1 nCoord nCP]
    
% Distance to polytope faces
    dp = pkg.geometry.distance.toPlane(P,O,N) ; % [nP 1 nCP]
    
% Approximation by faces ? (exact inside)
    if approx ; d = max(dp,[],3) ; return ; end

% Initialize
    d = NaN(size(P,1),1) ;
    
% INSIDE: use the "closest" face (maximum signed distance to face)
    isin = all(dp<=0,3) ; % [nP 1 nCP]
    d(isin) = max(dp(isin,:,:),[],3) ;
    
% OUTSIDE: sort cutting plane from the farthest to the closest
    dp = dp(~isin,:,:) ;
    dp = sort(dp,3,'descend') ;
    dp(dp<0) = NaN ; % do not take inside planes into account
    
% Take the closest feature
% 3D example: closest between
%   - the farthest plane
%   - the edge defined by the two farthest planes
%   - the vertice defined by the three farthest planes
    nMaxPl = min(nCP,nCoord) ;
%     dp = dp(:,:,1:nMaxPl) ;
%     dp2 = cumsum(dp.^2,2,'omitnan') ;
%     dp = sqrt(min(dp2,[],2)) ;

% Number of intersecting planes needed to form vertices
    dp = dp(:,:,1:nMaxPl) ;
    
    d(~isin) = sqrt(sum((dp.*(dp>0)).^2,3,'omitnan')) ;
    %d(~isin) = sum(dp>0,3) ;
end

function test
%% TEST THE DISTANCE
    nCoord = 2 ; nCP = 4 ; nP = 1000000 ;
    switch 1
        case 1 % random polytope
            O = 2*rand(nCP,nCoord)-1 ; 
            N = 2*rand(nCP,nCoord)-1 ;
        case 2 % centred random polytope
            O = 2*rand(nCP,nCoord)-1 ; 
            N = O ;
        case 3 % regular polytope
            O = pkg.math.combinations(repmat({[-1 1]},[1 nCoord])) ;
            N = O ;
    end
   
    N = N.*sign(sum(O.*N,2)) ;
    N = N./sqrt(sum(N.^2,2)) ;
    
    P = 4*(rand(nP,nCoord)-0.5)*2 ;
 %
    profile on
    tic ; d = pkg.geometry.distance.toPolytope(P,O,N,false) ; toc
    profile off
    d(d<0) = NaN ;
    
    cla ; axis equal ; axis tight
    P3D = padarray(P,[0 max(0,3-nCoord)],0,'post') ;
    s = scatter3(P3D(:,1),P3D(:,2),P3D(:,3)+d*0,5,d,'.') ;
    s.MarkerEdgeAlpha = .01 ;
    %s.MarkerEdgeAlpha = 'flat' ; s.AlphaData = double(abs(d)<=.025) ;
    colormap(jet(50))
    O3D = padarray(O,[0 max(0,3-nCoord)],0,'post') ;
    N3D = padarray(N,[0 max(0,3-nCoord)],0,'post') ;
    quiver3(O3D(:,1),O3D(:,2),O3D(:,3),N3D(:,1),N3D(:,2),N3D(:,3),'k') ;
    T3D = permute(O3D + 4*cat(3,-1,1).*N3D(:,[2 1 3]).*[-1 1 1],[3 1 2]) ;
    plot3(T3D(:,:,1),T3D(:,:,2),T3D(:,:,3),'-.k','linewidth',1)

end

