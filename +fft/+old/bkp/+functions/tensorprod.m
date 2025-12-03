function C = tensorprod(A,B,order)
% TENSORPROD of two gridded tensor field variables at a given order

    C = pkg.fft.old.Variable(A.Grid,'Size',[A.Size(1:end-order) B.Size(order+1:end)]) ;
    Ma = A.Size(1:end-order) ; 
    Na = A.Size(end-order+1:end) ; 
    Mb = B.Size(1:order) ;
    Nb = B.Size(order+1:end) ; 

    Ad = reshape(A.Data,[prod(Ma) prod(Na) 1 A.Grid.N]) ;
    Bd = reshape(B.Data,[prod(Mb) prod(Nb) 1 B.Grid.N]) ;

    Cd = pagemtimes(Ad,Bd) ; 

    C.Data = reshape(Cd,[C.Size C.Grid.N]) ; 

end

