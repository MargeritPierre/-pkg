function DATA = interp(DATA,IDX,ORDER,EXTRAP)
% Interpolation of a (N+P)D-array DATA with real-valued ND-indices IDX
% ORDER is the interpolation order <TODO>
% DATA [n1 n2 ... N p1 p2 ... P]
% IDX [nIdx Nt]
% data [nIdx p1 p2 ... P]
%
% The interpolation is achieved using a polynomial kernel:
%   Pn(x) = Xn(x).a (for any coordinates x)
%       with Xn(x) = [x^0 x^1 ... x^n]
%       and a = [a1 ; a2 ; ... ; an]
%   we know that Pn(idx) = Xn(idx).a = b
%   where b is the data values at integer indices idx
%   so a can be determined as a = Xn(idx)\b

% Default argumets
    if nargin<3 ; ORDER = 1 ; end
    if nargin<4 ; EXTRAP = true ; end

% Size information
    SZ = [size(DATA) 1 1] ;
    [nIdx,N] = size(IDX) ;
    n = SZ(1:N) ;
    p = SZ(N+1:end) ;
    
% kernel infos
    nk = ORDER+1 ;
    ord = 0:ORDER ;
    kk = (ord(:)-floor(ORDER/2)) ;

    
% Indices: integer part 'idx'
    if ORDER==0
        idx = round(IDX) ; % [nIdx N]
    else
        idx = floor(IDX) ; % [nIdx N]
    end

% Reshape the data array to a 2D matrix
    DATA = reshape(DATA,[prod(n) prod(p)]) ;

% Extrapolation
    imin = 1-min(kk) ;
    imax = n-max(kk) ;
    if EXTRAP 
        idx = max(imin,min(idx,imax)) ; % [nIdx N]
    else
        validIDX = ~any( IDX>imax | IDX<imin ,2) ;
        idx = idx(validIDX,:) ; % [nValidIdx N]
    end


% ZERO-TH ORDER: NEAREST NEIGHTBOR
    if ORDER==0
    % Linear indices
        idx = sum((idx-1).*[1 cumprod(n(1:end-1),2)],2) + 1 ; % [nValidIdx 1]
    % Data extraction
        validData = DATA(idx,:) ; % [nValidIdx prod(p)]
        
        
% HIGHER-ORDERS: USE THE INTERPOLATION KERNEL
    else
    % Extrapolation data
        if EXTRAP
            nValidIdx = nIdx ;
        else
            IDX = IDX(validIDX,:) ; % [nValidIdx N]
            nValidIdx = sum(validIDX) ;
        end

    % Indices: remaining residual 'dIDX'
        dIDX = IDX-idx ; % [nValidIdx N]

    % Polynomial kernel
        Xn0 = 1 ; Xn = 1 ;
        selDim = eye(N) ;
        for dim = 1:N
            shifts = kk.*selDim(dim,:) ; % [nk N]
            idx = repelem(idx,nk,1) + repmat(shifts,(nk^(dim-1))*nValidIdx,1) ; % [(nk^dim)*nValidIdx N]
            Xn0 = kron(Xn0,kk.^ord) ; % [nk^dim nk^dim]
            Xn = repelem(Xn,1,nk).*repmat(dIDX(:,dim).^ord,1,nk^(dim-1)) ; % [nValidIdx nk^dim]
        end

    % Linear indices
        idx = sum((idx-1).*[1 cumprod(n(1:end-1),2)],2) + 1 ; % [(nk^N)*nValidIdx 1]
        
    % vectors b
        idx = reshape(idx,[nk^N nValidIdx]) ; % [nk^N nValidIdx]
        b = reshape(DATA(idx,:),[nk^N nValidIdx*prod(p)]) ; % [nk^N nValidIdx*prod(p)]
    
    % amplitudes a
        a = Xn0\b ; % [nk^N nValidIdx*prod(p)]
        a = reshape(a,[nk^N nValidIdx prod(p)]) ;
        a = permute(a,[2 1 3]) ; % [nValidIdx nk^N prod(p)]
        
    % values
        validData = sum(Xn.*a,2) ; % matrix product Xn*a [nValidIdx 1 prod(p)]
        validData = reshape(validData,[nValidIdx prod(p)]) ;
        
    

    end


% Copy the data with invalid indices
    if EXTRAP
        DATA = validData ;
    else
        DATA = NaN([nIdx prod(p)]) ;
        DATA(validIDX,:) = validData ;
    end
    
% Reshape
    DATA = reshape(DATA,[nIdx p]) ;


end





function test
%% UNIT TEST, 1D

    t = 2*pi*linspace(0,1,10)' ;
    DATA = [sin(t) cos(2*t)] ;
    
    IDX = 1.1*linspace(0,1,1000)'*(numel(t)-1) + 1 ;
    ORDER = 5 ; 
    EXTRAP = true ;
    
    iDATA = pkg.data.interp(DATA,IDX,ORDER,EXTRAP) ;
    
    clf 
    plot(DATA,'o') ;
    plot(IDX,iDATA) ;
    
%% UNIT TEST 2D

    DATA = sin((1:10)'/1).*cos((1:12)/2) ;
    [iII,iJJ] = ndgrid(linspace(1,size(DATA,1),100),linspace(1,size(DATA,2),100)) ;
    IDX = [iII(:) iJJ(:)] ;
    ORDER = 4 ; 
    EXTRAP = true ;
   
    iDATA = pkg.data.interp(DATA,IDX,ORDER,EXTRAP) ;
    iDATA = reshape(iDATA,size(iII)) ;
    
    cla 
    surf(DATA,'Marker','o'...
            ,'MarkerEdgeColor','k'...
            ,'FaceColor','none'...
            ,'EdgeColor','none'...
            ) ;
    surf(iJJ,iII,iDATA,'facecolor','interp','edgecolor','none')
    %plot3(IDX(:,2),IDX(:,1),iDATA,'x')
    
    
%% SPEED TEST
    n = [1 1]*1000 ; N = numel(n) ;
    p = [1 1] ;
    nIdx = 10000000 ;
    ORDER = 0 ; 
    EXTRAP = true ;
    
    DATA = rand([n p]) ;
    IDX = rand(nIdx,N).*(n-1) + 1 ;
    
    profile on
    tic ; iDATA = pkg.data.interp(DATA,IDX,ORDER,EXTRAP) ; toc
    if N==2 ; tic ; iiData = interp2(DATA,IDX(:,2),IDX(:,1),'nearest') ; toc ; end
    profile off
    
    

end




