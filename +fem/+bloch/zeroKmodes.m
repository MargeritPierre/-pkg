function [U,omega] = zeroKmodes(K00,M,nModes)
% Compute the cut frequencies associated to zero wavenumber
% Solve (K00-omega^2*M)

    nDOFs = size(M,1) ;
    [U,w2] = eigs(K00,M,nModes,'sm') ;
    omega = sqrt(diag(w2)) ;
    


end