classdef sparse < matlab.mixin.CustomDisplay
%SPARSE nd-sparse array reprezentation
    
%% CLASS PROPERTIES
properties (Hidden)
    IDX % nd-indices of non-zero matrices
    VAL(:,1) % corresponding values
    SZ(1,:) = [0 0] % size of the array
end
   
%% CLASS CONSTRUCTOR(S)
methods
    function this = sparse(varargin)
    % Class Constructor
        switch nargin
            case 0
                return ;
            case 1 % Sparse(SZ), varargin sould contain the size only (full of zeros)
                this.SZ = floor(varargin{1}) ;
            case 2 % Sparse(IDX,VAL)
                this.IDX = floor(varargin{1}) ;
                this.VAL = varargin{2} ;
                this.SZ = max(this.IDX,[],1) ;
                this = accum(this) ;
            case 3 % Sparse(IDX,VAL,SZ)
                this.IDX = floor(varargin{1}) ;
                this.VAL = varargin{2} ;
                this.SZ = varargin{3} ;
                this = accum(this) ;
        end
    end
end

%% ACCUMULATION FUNCTION
methods
    function this = accum(this)
    % linear indices
        idx = sum((this.IDX-1).*cumprod([1 this.SZ(1:end-1)]),2) ;
    % sorting
        [idx,is] = sort(idx) ;
        val = this.VAL(is) ;
    % unique
        isNew = [true ; ~logical(diff(idx))] ;
        this.IDX = this.IDX(isnew,:) ;
    % summation
        val = cumsum(val) ;
        val = val(cirshift(isNew,-1)) ;
        this.VAL = [val(1) ; diff(val)] ;
    end
end

%% OVERRIDE BUILTIN ARRAY INFO FUNCTIONS
methods
    function sz = size(this,dim)
        if nargin<2 ; dim = 1:numel(this.SZ) ; end
        if dim>numel(this.SZ) ; sz = 1 ; return ; end
        sz = this.SZ(dim(:)') ;
    end
    function n = numel(this) ; n = prod(this.SZ) ; end
    function n = nnz(this) ; n = numel(this.VAL) ; end
end

%% SUBSREF,SUBSASSIGN
methods
    function this = subsref(this,S)
        
    end
    
    function this = subsasgn(this,S,B) 
    end
end

%% OBJECT DISPLAY FUNCTIONS
methods (Access = protected)
    function header = getHeader(this)
        className = matlab.mixin.CustomDisplay.getClassNameForHeader(this) ;
        sizeStr = arrayfun(@num2str,this.SZ,'UniformOutput',false) ;
        header = [strjoin(sizeStr,'x') ' ' className ' ' class(this.VAL)] ;
        header = sprintf('\t%s\n',header) ;
    end
end


%% END ===============
end

