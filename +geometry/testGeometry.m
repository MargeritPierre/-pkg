%% UNIT TESTS for the GEOMETRY package

clf
clear all
clc

%set(gca,'toolbar',[])
set(gca,'interactions',[]) ;
set(gca,'xminortick','on','xminorgrid','on','yminortick','on','yminorgrid','on','zminortick','on','zminorgrid','on')
set(gca,'boxstyle','full','box','on')
grid on

    % Creation of base geometries
        p = pkg.geometry.Point ;
        p.Position = [0 0 0] ;
        p.draw ;
        
        
        
        