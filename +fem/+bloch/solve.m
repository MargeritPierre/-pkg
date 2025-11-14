%% BLOCH WAVE COMPUTATION
% (K00 + mu.K1 + mu^2.K2 - w^2.M).U = 0
% V = mu.U
% ([w^2.M-K0 0;0 K2] - mu.[K1 K2;K2 0])[U;V] = [0;0]

function [K,U,omega] = solve(K00,K0i,Kij,M,freq,dir,nModes)

    omega = 2*pi*freq ;
    dir(isnan(dir)) = 0 ; % avoid multiplication by NaNs
    nDir = size(dir,1) ;
    nFreq = numel(omega) ;
    nDOFs = size(M,1) ;

    vK0i = cellfun(@(M)reshape(M,[],1),K0i,'uni',false) ;
    vK0i = cat(2,vK0i{:}) ;
    vKij = cellfun(@(M)reshape(M,[],1),Kij,'uni',false) ;
    vKij = cat(2,vKij{:}) ;
    I = speye(nDOFs) ; O = sparse(nDOFs,nDOFs) ;
    if ~issparse(M) ; I = full(I) ; O = full(O) ; end

    if nargin<7 || isempty(nModes) ; nModes = 2*nDOFs ; end
    K = NaN(nModes,nFreq,nDir) ;
    U = NaN(nDOFs,nModes,nFreq,nDir) ;
    
    wtbr = waitbar(0,'Bloch Wave Computation..') ;
    for dd = 1:nDir
    % Projected stiffnesses
        n = dir(dd,:) ;
        n(end+1:numel(K0i)) = 0 ;
        nn = n(:).*n(:)' ; 
        K1 = reshape(vK0i*sparse(n(:)),[nDOFs nDOFs]) ;
        K2 = reshape(vKij*sparse(nn(:)),[nDOFs nDOFs]) ;
    % For each frequency..
        for ww = 1:nFreq
            K0 = omega(ww)^2*M-K00 ;
            A = blkdiag(K0,I) ; 
            B = [K1 K2 ; I O] ;
            if nModes==2*nDOFs % extract all modes
                [uu,ku] = eig(A,B) ;
            else
                [uu,ku] = eigs(A,B,nModes,'sm') ;
            end
            K(:,ww,dd) = diag(ku) ;
            U(:,:,ww,dd) = reshape(uu(1:end/2,:),nDOFs,nModes) ;
            wtbr = waitbar((ww+nFreq*(dd-1))/nFreq/nDir,wtbr) ;
        end
    end
    delete(wtbr) ;

end