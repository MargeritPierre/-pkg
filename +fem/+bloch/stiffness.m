%% MATERIAL STIFFNESS MATRIX
function C = stiffness(E,G)
    E = reshape(E,1,1,[]) ; G = reshape(G,1,1,[]) ; O = zeros(size(E)) ;
    nu = E./(2*G)-1 ;
    C = E./(1+nu)./(1-2*nu).*[...
                              1-nu nu nu O O O ; ...
                              nu 1-nu nu O O O ; ...
                              nu nu 1-nu O O O ; ...
                              O O O .5-nu O O ; ...
                              O O O O .5-nu O ; ...
                              O O O O O .5-nu ; ...
                             ] ;
end