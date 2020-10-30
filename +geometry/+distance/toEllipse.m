function d = toEllipse(P,CENTER,SEMIAXES,ROTATION)
%TOELLIPSE SIGNED Distance from a list of points P to an Ellipse
% see https://wet-robots.ghost.io/simple-method-for-distance-to-ellipse/
% inputs:
%   - P [nP nCoord]
%   - CENTER [1 nCoord 1] center of the ellipse
%   - SEMIAXES [1 nCoord 1] semiaxes of the ellipse
%   - ROTATION [1 nCoord 1] rotation of the ellipse (radians)
% output: 
%   - d [nP 1 1] signed distance (d<0 inside, d>0 outside)

    if isempty(P) ; d = [] ; return ; end
    
    % We start by centering the ellipse on origin
        P = P-CENTER ;
    % Then we rotate the frame of -rotation
        P = P*[cos(ROTATION) -sin(ROTATION) ; sin(ROTATION) cos(ROTATION)] ;
    % Any point pe on the ellipse is then :
        a = SEMIAXES(1) ; b = SEMIAXES(2) ;
        pe = @(t)[a*cos(t) b*sin(t)] ;
    % A point on the "evolute" is:
        pv = @(t)(b^2-a^2).*[-cos(t).^3./a sin(t).^3./b] ;
    % We push every point to the first quadrant
        P = abs(P) ;
    % And the initial guess is taken as follows:
        t = pi/4 ; %atan2(a*P(:,2),b*P(:,1)) ;
    % Angle tolerance
        tol = 1e-6 ;
    % Maximum number of iterations
        maxIt = 5 ;
    % Newton descent
        dt = Inf ;
        it = 0 ;
        while it<maxIt && max(abs(dt))>tol
            it = it+1 ;
            %
            x = pe(t) ;
            e = pv(t) ;
            %
            r = x - e ;
            q = P - e ;
            %
            R = sqrt(sum(r.^2,2)) ;
            Q = sqrt(sum(q.^2,2)) ;
            %
            dc = R.*asin( (r(:,1).*q(:,2) - r(:,2).*q(:,1)) ./ (R.*Q) ) ;
            dt = dc./sqrt(a^2 + b^2 - sum(x.^2,2)) ;
            %
            t = t + dt ;
            t = min(pi/2,max(0,t)) ;
        end
    % Restablish to the four quadrants
        x = pe(t) ;
    % Absolute distance
        d = sqrt(sum((x-P).^2,2)) ;
    % Signed Distance
        d = d.*sign(P(:,1).^2/a^2+P(:,2).^2/b^2-1) ; 

end

