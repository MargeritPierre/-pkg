function [p,d] = closestToMouse(pts,ax)
% Get the closest data point to the mouse pointer
    if nargin<2 ; ax = gca ; end
% Mouse position in axes 
    line = ax.CurrentPoint ;
% Set to limits if rotation perpendicular
    lims = [ax.XLim(:) ax.YLim(:) ax.ZLim(:)] ;
    line(isinf(line)) = lims(isinf(line)) ;
% Set scales if needed
    islog = ismember({ax.XScale ax.YScale ax.ZScale},'log') ;
    line(:,islog) = log10(line(:,islog)) ;
    lims(:,islog) = log10(lims(:,islog)) ;
    pts(:,islog,:) = log10(pts(:,islog,:)) ;
% Normalize with limits
    m = min(lims,[],1) ;
    R = range(lims,1) ;
    line = (line-m)./R ;
    pts = (pts-m)./R ;
% Distance from line to all points
    u = diff(line,1,1)./sqrt(sum(diff(line,1,1).^2,2)) ;
    v = pts-line(1,:) ;
    t = sum(v.*u,2) ;
    d = sqrt(abs(sum(v.^2,2) - t.^2)) ;
    [d,p] = min(d(:))  ;
end