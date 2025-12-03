clear all

N = [(2^6)*[2 2 1] 3 3 ] ;

fcn = @(A,dims)@()fftd(A,dims) ;
% fcn = @(A,dims)@()pagemtimes(A,'none',A,'ctranspose') ;
% fcn = @(A,dims)@()sum(A,dims) ;
% fcn = @(A,dims)@()sum(abs(A).^2,dims) ;
% fcn = @(A,dims)@()A+A ;
% fcn = @(A,dims)@()arrayfun(@plus,A,A) ;
%fcn = @(A,dims)@()permute(A,[dims setdiff(1:ndims(A),dims)]) ;

dims = [...
            ...num2cell(1:numel(N),1)' ; ...
            {1:numel(N)} ; ...
            {find(N==max(N))} ; ... % gridfun
            {find(N==min(N))} ; ... % tensorfun
        ] ;

A = single(randn(N) + 1i*randn(N)) ;
Agpu = gpuArray(complex(A)) ;

T = table ;
T.Dim = dims ;
T.N = cellfun(@(d)N(d),dims,'uni',false) ;

disp([])
disp("N="+mat2str(N))
for d = 1:numel(dims)
    dd = dims{d} ;
    % CPU
    T.CPU(d) = timeit(fcn(A,dd),1) ;
    disp("CPU/"+mat2str(dd)+":"+string(T.CPU(d)))
    % GPU
    T.GPU(d) = gputimeit(fcn(Agpu,dd),1) ;
    disp("GPU/"+mat2str(dd)+":"+string(T.GPU(d)))
    % Check error
    tic ; cpuFcn = fcn(A,dd) ; %toc
    tic ; gpuFcn = fcn(Agpu,dd) ; %toc
    T.RelError(d) = norm(reshape(cpuFcn()-gather(gpuFcn()),[],1))/norm(reshape(cpuFcn(),[],1)) ;
end

T.Speedup = T.CPU./T.GPU ;
disp(T)


function A = fftd(A,dims)
    % if isa(A,'gpuArray') 
    %     perm = [dims setdiff(1:ndims(A),dims)] ; 
    %     % perm = [setdiff(1:ndims(A),dims) dims] ; 
    %     A = permute(A,perm) ;
    %     dims = 1:numel(dims) ;
    % end
    for d = dims
        A = fft(A,[],d) ;
        A = ifft(A,[],d) ;
    end
    % if isa(A,'gpuArray') 
    %     [~,is] = sort(perm) ;
    %     A = permute(A,is) ; 
    % end
end