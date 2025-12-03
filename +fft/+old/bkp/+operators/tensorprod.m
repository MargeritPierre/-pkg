classdef tensorprod < pkg.fft.old.operator
% TENSORPROD Defines the tensor product between two gridded tensor fields

    properties
        Order(1,1) uint32
    end
    
    methods
        function this = tensorprod(A,B,order)
        % Constructor
            this.Inputs = [A B] ;
            this.Order = order ;
            this.Outputs = pkg.fft.old.field(A.Grid,'Size',[A.Size(1:end-order) B.Size(order+1:end)]) ;
        end

        function apply(this)
            A = this.Inputs(1) ; B = this.Inputs(2) ; C = this.Outputs ;
            nA = A.Size(1:end-this.Order) ; nB = B.Size(this.Order+1:end) ; 
            Ad = reshape(A.Data,prod(nA),[],prod(C.Grid.N)) ;
            Bd = reshape(B.Data,[],prod(nB),prod(C.Grid.N)) ;
            Cd = pagemtimes(Ad,Bd) ; 
            Cd = reshape(Cd,[C.Size C.Grid.N]) ; 
            C.Data = Cd ;
        end
    end
end

