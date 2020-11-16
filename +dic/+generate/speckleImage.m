function I = speckleImage(MASK,nPts,medFiltSize,gaussFiltSize,edgeWidth,backgroundLvl)
%SPECKLEIMAGE Generate a speckle image from a binary mask representing the object
%geometry
% MASK: binary mask [nI nJ]
% nPts: number of initial speckle points
% medFiltSize: median filter radius
% gaussFiltSize: gaussian filter radius
% edgeWidth: add white edges with a given width
% backgroundLvl: level of the background

% Default mask (for demo)
if nargin<1 
    sz = 512*[1 1] ; % !!! [nJ nI] !!!
    [JJ,II] = meshgrid(1:sz(1),1:sz(2)) ;
    shape = ... [0 0 ; 1 0 ; 1 1 ; 0 1].*sz.*[8 6]/10 + [1 2]/10.*sz ... rectangle
            real(exp(2i*pi*(0:1/7:1-1e-9)').*[1 -1i]).*sz.*[4 4]/10 + sz/2 ... polyhedron
            ;
    MASK = inpolygon(JJ(:),II(:),shape(:,1),shape(:,2)) ;
    MASK = reshape(MASK,flip(sz)) ; % [nI nJ]    
end
sz = size(MASK) ;

% Other default inputs
if nargin<2 ; nPts = round(7/10*prod(sz)) ; end
if nargin<3 ; medFiltSize = 7 ; end
if nargin<4 ; gaussFiltSize = 9 ; end
if nargin<4 ; edgeWidth = 3 ; end
if nargin<4 ; backgroundLvl = 0.1 ; end

% Texture generation
I = false(sz) ;
% Random distribution of active pixels
    ii = randi([1 sz(1)],nPts,1) ;
    jj = randi([1 sz(2)],nPts,1) ;
    I(sub2ind(size(I),ii,jj)) = true ;
% Median Filtering
    I = medfilt2(I,[1 1]*medFiltSize) ;
            
% Apply the shape
I = double(I) ;
% Set the background
    I(~MASK) = backgroundLvl ;
% Add edges
    I(MASK & ~bwmorph(MASK,'thin',edgeWidth)) = 1 ;

% Gaussian filter
% Kernel
    f = gausswin(gaussFiltSize) ; 
    f = f(:)*f(:)' ; 
    f = f./sum(abs(f(:))) ;
% Convolution
    I = conv2(I,f,'same') ;
    
end

%% UNIT TESTS
function test
%% GENERATE DEMO
    
    I = pkg.dic.generate.speckleImage ;

    % Display
        clf
        imagesc(repmat(I,[1 1 3]))
        axis tight, axis equal, axis ij, box on
end



