function var = fft(var)
% FFT of a gridded tensor field variable
    if var.Grid.Fourier ; return ; end
    if nargout>0 ; var = copy(var) ; end % might do in-place computations if possible
    var.Grid = fft(var.Grid) ;
    if ~isempty(var.Data)
        for d = 1:var.Grid.ndims
            var.Data = fft(var.Data,var.Grid.N(d),var.Order+d) ;
        end
    end
end

