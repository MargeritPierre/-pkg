%% DEMO SCRIPT FOR VIEWPORTS

%% CREATE A DEFAULT VIEWPORT

clearvars % delete previously defined viewports
close all

vp = pkg.ui.viewport.Viewport ;
vp.Interactions = 'all' ;

P = rand(1000,3) ;
plot3(P(:,1),P(:,2),P(:,3),'.k','markersize',10)
%axis equal
xlabel x
ylabel y
zlabel z