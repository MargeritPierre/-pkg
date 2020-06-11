function R = rotmat(ANGLE,AXIS)
%ROTATIONMATRIX return the rotation matrix corresponding to the rotation of
%an ANGLE (in radians) about an AXIS
% AXIS: [nRot nCoord==3] ;
% ANGLE = [nRot 1] ;
% see https://fr.wikipedia.org/wiki/Matrice_de_rotation
if nargin<2 ; AXIS = [0 0 1] ; end

AXIS = AXIS(:,:) ;
ANGLE = ANGLE(:) ;

[nRot,nCoord] = size(AXIS(:,:)) ;

if nCoord~=3 ; error('This function wors in 3D only') ; end
if nRot==1 ; nRot = numel(ANGLE) ; end
if numel(ANGLE)~=1 && nRot~=numel(ANGLE)
    error('Cannot process the given inputs') ;
end

AXIS = AXIS./sqrt(sum(AXIS.^2,2)) ;

AXIS = permute(AXIS,[3 2 1]) ;
ANGLE = permute(ANGLE,[3 2 1]) ;

P = AXIS.*permute(AXIS,[2 1 3]) ;
I = eye(3) ;
O = AXIS(1,1,:).*0 ;
Q = [   O -AXIS(:,3,:) AXIS(:,2,:) ; ...
        AXIS(:,3,:) O -AXIS(:,1,:) ; ...
        -AXIS(:,2,:) AXIS(:,1,:) O ; ] ;
    
R = P + cos(ANGLE).*(I-P) + sin(ANGLE).*Q ;

end

