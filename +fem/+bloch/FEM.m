%% BLOCH WAVE PROPAGATION SYSTEM MATRICES
% u(x,t) = U(x)*exp(1i*omega*t - 1i*k*x)
% (K00 + k.K1 + k^2.K2 - w^2.M).U = 0

function [K00,K0i,Kij,M,P] = FEM(mesh,C,rho,sig,perVec,nUCoord,nKCoord)

    if nargin<6 ; nUCoord = 3 ; end
    if nargin<7 ; nKCoord = 3 ; end

    % Mesh interpolation
    [ee,we,ie] = mesh.integration() ;
    nQP = numel(we) ;
    W = diag(sparse(we)) ;
    O = sparse(nQP,mesh.nNodes) ;
    N = mesh.interpMat(ee,ie) ;
    D = mesh.diffMat(ee,ie) ; 
    [D{end+1:3}] = deal(O) ;

    % EPS = [E11;E22;E33;2E13;2E23;2E12] = B*U = (B0 + sum(B{i}*k_i))*u = [
    B0 = [D{1} O O ; O D{2} O ; O O D{3} ; D{3} O D{1} ; O D{3} D{2} ; D{2} D{1} O] ;
    Bi = cell(3,1) ;
    Bi{1} = -1i*[N O O ; O O O ; O O O ; O O N ; O O O ; O N O] ;
    Bi{2} = -1i*[O O O ; O N O ; O O O ; O O O ; O O N ; N O O] ;
    Bi{3} = -1i*[O O O ; O O O ; O O N ; N O O ; O N O ; O O O] ;

    % DISTRIBUTE MATERIAL STIFFNESS & WEIGHTS
    if size(C,3)==1 % homogeneous solid
        CW = kron(sparse(C),W) ;
    else % heterogeneous
    % material stiffness at quadrature points
        if size(C,3)==nQP
            Cqp = C ;
        elseif size(C,3)==mesh.nElems
            Cqp = C(:,:,ie) ;
        elseif size(C,3)==mesh.nNodes
            Cqp = reshape((N*reshape(C,[6*6 mesh.nNodes]).').',[6 6 nQP]) ;
        end
        Cqp = Cqp.*reshape(we,1,1,[]) ; % times weights
        iii = (0:5)'*nQP + zeros(1,6) + reshape(1:nQP,1,1,[]) ;
        jjj = zeros(6,1) + (0:5)*nQP + reshape(1:nQP,1,1,[]) ;
        CW = sparse(iii(:),jjj(:),Cqp(:),6*nQP,6*nQP) ;
    end

    % STIFFNESS MATRICES
    K00 = B0'*CW*B0 ; % zero-order matrix
    K0i =  cell(nKCoord,1) ; % first-order matrices
    Kij = cell(nKCoord,nKCoord) ; % second order matrices
    for ii = 1:nKCoord 
        K0i{ii} = B0'*CW*Bi{ii} + Bi{ii}'*CW*B0 ;
        for jj = 1:nKCoord
            Kij{ii,jj} = Bi{ii}'*CW*Bi{jj} ; % Kij{ii,jj} = Kij{jj,ii}'
        end
    end

    % MASS MATRICES
    if numel(rho)==1
        rhoW = rho*we ;
    elseif numel(rho)==nQP
        rhoW = rho(:).*we ;
    elseif numel(rho)==mesh.nNodes
        rhoW = (N*rho(:)).*we ;
    elseif numel(rho)==mesh.nElems
        rhoW = rho(ie).*we ;
    end
    rhoW = diag(sparse(rhoW)) ;
    M = (N'*rhoW*N) ;
    M = blkdiag(M,M,M) ;
    
    % PRESTRESS INFLUENCE
    if nargin>=4 && any(abs(sig)>eps)
    % Stress components on quadrature points
        if size(sig,1)==1 % homogeneous stress
            sigW = we.*sig ;
        elseif size(sig,1)==mesh.nNodes % nodal stress
            sigW = W*(N*sig) ;
        else % element stress
            sigW = W*sig(ie,:) ;
        end
    % Distibute as a sparse matrix
        sigW = permute(sigW(:,[1 6 4;6 2 5;4 5 3]),[2 3 1]) ; % (3,3,nQP) 
        iii = (0:2)'*nQP + zeros(1,3) + reshape(1:nQP,1,1,[]) ;
        jjj = zeros(3,1) + (0:2)*nQP + reshape(1:nQP,1,1,[]) ;
        sigW = sparse(iii(:),jjj(:),sigW(:),3*nQP,3*nQP) ;
    % Displacement gradients [u1,1 ; u1,2 ; u1,3 ; ... ; u3,3]
        G0 = cat(1,D{:}) ;
        Gk = -1i*N ;
    % Update stiffness matrices
        k00 = G0'*sigW*G0 ;
        K00 = K00 + blkdiag(k00,k00,k00) ;
        for ii = 1:nKCoord 
            Gi = kron(sparse(ii,1,1,3,1),Gk) ;
            k0i = G0'*sigW*Gi ;
            k0i = k0i + k0i' ; % symmetric hermitian
            K0i{ii} = K0i{ii} + blkdiag(k0i,k0i,k0i) ;
            for jj = 1:nKCoord
                Gj = kron(sparse(jj,1,1,3,1),Gk) ;
                kij = Gi'*sigW*Gj ;
                kij = kij + kij' ; % symmetric hermitian
                Kij{ii,jj} = Kij{ii,jj} + blkdiag(kij,kij,kij) ;
            end
        end
    end

    % PERIODICITY MATRICES
    if nargin<5 ; perVec = diag(range(mesh.boundingBox,1)) ; end
    if isempty(perVec) ; P = speye(mesh.nNodes) ;
    else ; P = mesh.perNodeMat(perVec) ;
    end
    % THEN ...
    P(:,sum(P,1)==0) = [] ; % delete unused DOFs
    P = repmat({P},[1 nUCoord]) ;
    P = blkdiag(P{:}) ; % periodicity on the 3 coordinates (constant coordinate frame)
    P3D = P ; if nUCoord<3 ; P3D(end+1:3*mesh.nNodes,:) = 0 ; end
    
    % APPLY PERIODICITY
    M = P3D'*M*P3D ;
    K00 = P3D'*K00*P3D ;
    for ii = 1:nKCoord
        K0i{ii} = P3D'*K0i{ii}*P3D ;
        for jj = 1:nKCoord
            Kij{ii,jj} = P3D'*Kij{ii,jj}*P3D ;
        end
    end
    
end