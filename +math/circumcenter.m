function [P,R] = circumcenter(TRI)
%CIRCUMCENTER Return the circumcenters of a collection of triangles
% TRI: [3 nCoord nTris]
% P: [nTris nCoord] ;
% R: [nTris 1]
% see https://en.wikipedia.org/wiki/Circumscribed_circle

sz = size(TRI) ;
if sz(1)~=3 ; error('The input TRI should have the size [3 nCoord nTris]') ; end

% Place A at the origin
    C = TRI(end,:,:) ;
    TRI = TRI-C ;
    
% Circumcenter coordinates
    norm2A = sum(TRI(1,:,:).^2,2) ;
    norm2B = sum(TRI(2,:,:).^2,2) ;
    AscalB = sum(TRI(1,:,:).*TRI(2,:,:),2) ;
    norm2AcrossB = norm2A.*norm2B - AscalB.^2 ;
    E = norm2A.*TRI(2,:,:) - norm2B.*TRI(1,:,:) ;
    EscalA = sum(E.*TRI(1,:,:),2) ;
    EscalB = sum(E.*TRI(2,:,:),2) ;
    P = (EscalB.*TRI(1,:,:) - EscalA.*TRI(2,:,:))./(2*norm2AcrossB) ;
    
% Reshape
    P = permute(P,[3 2 1]) ;
    R = sqrt(sum(P.^2,2)) ;
    
% Add Origin back
    P = P + permute(C,[3 2 1]) ;

end

