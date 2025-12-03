function [outGPU,outHost] = gpuBench()
%GPUBENCH  MATLAB GPU Benchmark
%   GPUBENCH times different MATLAB GPU tasks and compares the execution
%   speed with the speed of several other GPUs.  The tasks are:
%
%    Backslash   Matrix left-division.    Floating point, regular memory access.
%    MTimes      Matrix multiplication.   Floating point, regular memory access.
%    FFT         Fast Fourier Transform.  Floating point, irregular memory access.
%
%   Each task is run for a range of array sizes and the results are tabulated
%   in an HTML report.  GPUBENCH can take several minutes to complete - please
%   be patient! Note that if your GPU is also driving your monitor then
%   the display may become unresponsive during testing.
%
%   GPUBENCH runs each of the tasks and shows a report indicating how the
%   current GPU compares to other systems.
%
%   T = GPUBENCH returns a data structure containing all of the results and
%   does not generate the report.
%
%   Fluctuations of up to ten percent in the measured times of repeated
%   runs on a single machine are not uncommon.  Your own mileage may vary.
%
%   This benchmark is intended to compare performance different GPUs on one
%   particular version of MATLAB.  It does not offer direct comparisons
%   between different versions of MATLAB.
%
%   See also: BENCH, gpuBenchReport

% Unused tasks:
%    Mandelbrot  Calculate a Mandelbrot Set.  Floating point, regular memory access.

%   Copyright 2011-2025 The MathWorks, Inc.

% These should probably be options in future.
[displayProgess,graphicalProgress] = gpubench.getDisplayMode();

% Check for the right MATLAB version and availability of PCT
gpubench.checkMATLABVersion();
gpubench.checkPCT();

% Check for a GPU. We give the option of running without a GPU so that
% users can evaluate what benefits a GPU might give.
hasGPU = parallel.gpu.GPUDevice.isAvailable();

% GPU benchmarking explicitly asked for. Error if GPU not available.
if nargout == 1 && ~hasGPU
    iThrowNoGPUError();
end

% Both CPU and GPU bench requested.
doHostOnNonGPUSystem = false;
if ~hasGPU
    title = 'Continue without a GPU?';
    question = ['The GPU could not be used. ' ...
        'Do you wish to continue and collect results for your CPU?'];
    buttons = {'Collect CPU results', 'Stop'};
    answer = questdlg( question, title, buttons{:}, buttons{end} );
    if ~strcmp(answer,buttons{1})
        iThrowNoGPUError();
    end
    doHostOnNonGPUSystem = true;
end

% Show the progress bar early as it is the first sign we are doing
% something!
% Do we need to measure the host stuff?
doHost = (nargout~=1) || doHostOnNonGPUSystem;
totalTasks = 6*(hasGPU + doHost);
title = sprintf('Running GPUBench v%s',gpubench.version);

if displayProgess
    fprintf('%s:\n',title);
end
if graphicalProgress
    wb = waitbar(0, title, 'Name', title);
    cleaner = onCleanup(@() close(wb));
else
    wb = [];
end

% Initialize the data object
release = regexp( version, 'R\d*[ab]', 'match' );
timestamp = now(); %#ok<TNOW1> For backwards compatibility keep as datenum
gpuData = gpubench.PerformanceData( ...
    release{1}, ...
    gpubench.cpuinfo(), ...
    gpubench.gpuinfo(), ...
    false, ... % isHostData
    timestamp );
hostData = gpubench.PerformanceData( ...
    release{1}, ...
    gpubench.cpuinfo(), ...
    struct(), ...
    true, ... % isHostData
    timestamp );

tasksDone = 0;
if hasGPU
    gpuData = runBackslash( gpuData, 'single', 'GPU', displayProgess );
    tasksDone = iMaybeUpdateWaitbar(graphicalProgress, tasksDone, totalTasks, wb);
    
    gpuData = runBackslash( gpuData, 'double', 'GPU', displayProgess );
    tasksDone = iMaybeUpdateWaitbar(graphicalProgress, tasksDone, totalTasks, wb);

    gpuData = runMTimes( gpuData, 'single', 'GPU', displayProgess );
    tasksDone = iMaybeUpdateWaitbar(graphicalProgress, tasksDone, totalTasks, wb);
    
    gpuData = runMTimes( gpuData, 'double', 'GPU', displayProgess );
    tasksDone = iMaybeUpdateWaitbar(graphicalProgress, tasksDone, totalTasks, wb);

    gpuData = runFFT( gpuData, 'single', 'GPU', displayProgess );
    tasksDone = iMaybeUpdateWaitbar(graphicalProgress, tasksDone, totalTasks, wb);

    gpuData = runFFT( gpuData, 'double', 'GPU', displayProgess );
    tasksDone = iMaybeUpdateWaitbar(graphicalProgress, tasksDone, totalTasks, wb);
end

if doHost
    hostData = runBackslash( hostData, 'single', 'CPU', displayProgess );
    tasksDone = iMaybeUpdateWaitbar(graphicalProgress, tasksDone, totalTasks, wb);

    hostData = runBackslash( hostData, 'double', 'CPU', displayProgess );
    tasksDone = iMaybeUpdateWaitbar(graphicalProgress, tasksDone, totalTasks, wb);

    hostData = runMTimes( hostData, 'single', 'CPU', displayProgess );
    tasksDone = iMaybeUpdateWaitbar(graphicalProgress, tasksDone, totalTasks, wb);

    hostData = runMTimes( hostData, 'double', 'CPU', displayProgess );
    tasksDone = iMaybeUpdateWaitbar(graphicalProgress, tasksDone, totalTasks, wb);

    hostData = runFFT( hostData, 'single', 'CPU', displayProgess );
    tasksDone = iMaybeUpdateWaitbar(graphicalProgress, tasksDone, totalTasks, wb);

    hostData = runFFT( hostData, 'double', 'CPU', displayProgess );
    iMaybeUpdateWaitbar(graphicalProgress, tasksDone, totalTasks, wb);
end

if nargout
    % User requested raw data
    outGPU = gpuData;
    outHost = hostData;
else
    % Produce report
    reportData = {};
    if hasGPU
        reportData{end+1} = gpuData;
    end
    if doHost
        reportData{end+1} = hostData;
    end
    % Make sure the progress bar closes before the report opens its own
    if graphicalProgress
        clear('cleaner')
    end
    web( gpuBenchReport( reportData{:} ) );
end


%-------------------------------------------------------------------------%
function data = runFFT( data, type, device, displayProgess )
% Work out the maximum size we should run
safetyFactor = 10; % Based on trial and error. Requiring 10x the input seems safe.
sizes = getTestSizes( type, safetyFactor, device );
times = inf( size( sizes ) );
avgTime = 0;

% Create a RandStream on CPU or GPU for getting data.
rs = iGetRandstream( 123, device );

if displayProgess
    fprintf( '* FFT (%s, %s)       ', device, type );
    cleanup = onCleanup( @() fprintf( ' done.\n' ) );
end

for ii=1:numel(sizes)
    % Check for getting close to time-out
    if tooCloseToTimeout( avgTime, device )
        times(ii) = nan;
        continue;
    end

    N = sizes(ii);
    try
        A = complex( rand( rs, N, 1, type ), rand( rs, N, 1, type ) );

        times(ii) = iTimeit( device, @()fft( A ) );
        avgTime = times(ii);

        if displayProgess
            fprintf('.');
        end

    catch err
    end
end

% Clear any dud results
sizes(isnan( times )) = [];
times(isnan( times )) = [];

data = addResult( data, 'FFT', type, sizes, 5*sizes.*log2(sizes), times );

%-------------------------------------------------------------------------%
function data = runMTimes( data, type, device, displayProgess )
safetyFactor = 3.5; % Space for two inputs plus one output and a bit to spare
sizes = getTestSizes( type, safetyFactor, device );

times = inf( size( sizes ) );
avgTime = 0;

% Create a RandStream on CPU or GPU for getting data.
rs = iGetRandstream( 123, device );

if displayProgess
    fprintf( '* Mtimes (%s, %s)    ', device, type );
    cleanup = onCleanup( @() fprintf( ' done.\n' ) );
end

N = round( sqrt( sizes ) );
for ii=1:numel(sizes)
    % Check for getting close to time-out
    if tooCloseToTimeout( avgTime, device )
        times(ii) = nan;
        continue;
    end
    try
        A = rand( rs, N(ii), N(ii), type );
        B = rand( rs, N(ii), N(ii), type );

        times(ii) = iTimeit( device, @()A*B );
        avgTime = times(ii);

        if displayProgess
            fprintf('.');
        end
    catch err
    end
end

% Clear any dud results
N(isnan( times )) = [];
times(isnan( times )) = [];

data = addResult( data, 'MTimes', type, N.*N, N.*N.*(2.*N-1), times );

%-------------------------------------------------------------------------%
function data = runBackslash( data, type, device, displayProgess )
safetyFactor = 1.5; % One full-sized matrix plus two vectors, so 1.5 is plenty
sizes = getTestSizes( type, safetyFactor, device );

% Limit the sizes to 1e8 for now to prevent problems
sizes(sizes>1e8) = [];
times = inf( size( sizes ) );
avgTime = 0;

% Create a RandStream on CPU or GPU for getting data.
rs = iGetRandstream( 123, device );

if displayProgess
    fprintf( '* Backslash (%s, %s) ', device, type );
    cleanup = onCleanup( @() fprintf( ' done.\n' ) );
end

N = round( sqrt( sizes ) );
for ii=1:numel( sizes )
    % Check for getting close to time-out
    if tooCloseToTimeout( avgTime, device )
        times(ii) = nan;
        continue;
    end

    try
        A = 100*eye( N(ii), N(ii), type ) + rand( rs, N(ii), N(ii), type );
        b = rand( rs, N(ii), 1, type );
        times(ii) = iTimeit( device, @()A\b );
        avgTime = times(ii);

        if displayProgess
            fprintf('.');
        end
    catch err
    end
end

% Clear any dud results
N(isnan( times )) = [];
times(isnan( times )) = [];

data = addResult( data, 'Backslash', type, N.*N, round( 2/3*N.^3 + 3/2*N.^2 ), times );

%-------------------------------------------------------------------------%
function sizes = getTestSizes( type, safetyFactor, device )
% Return the maximum number of elements that will fit in the device memory
elementSize = gpubench.sizeof( type );
% If no GPU to get memory size, so just go for 4GB
freeMem = 4*2^30;

if strcmpi( device, 'CPU' )
    % On the host everything takes longer, so don't go as far
    safetyFactor = safetyFactor*2;
else
    % Use as much memory as we can.
    gpu = gpuDevice();
    freeMem = gpu.FreeMemory;
end

maxNumElements = floor( freeMem / (elementSize*safetyFactor) );
if isnan( maxNumElements ) || maxNumElements < 1e6
    error( 'gpuBench:NotEnoughMemory', 'Not enough free device memory to run tasks' );
end

% We want powers of two up to this size
maxPower = floor( log2( maxNumElements ) );
sizes = power( 2, 10:2:maxPower );

%-------------------------------------------------------------------------%
function stopNow = tooCloseToTimeout( time, device )
% Should a test stop early to avoid triggering the device time-out?
stopNow = false;
if strcmpi( device, 'CPU' )
    % On the host there is no time limit
else
    gpu = gpuDevice();
    % If the kernel has a timeout it is typically 2-5 seconds. If we have
    % just done a size that takes (on average) more than a tenth of a second,
    % the next size will likely trigger the timeout.
    stopNow = (gpu.KernelExecutionTimeout && time>0.1);
end

%-------------------------------------------------------------------------%
function rs = iGetRandstream( seed, device )
% Get a Randstream that will generate data on CPU or GPU as appropriate
if verLessThan('MATLAB','9.6')
    % Prior to R2019a (v9.6) ThreeFry was not available so use CombRecursive.
    args = {'CombRecursive', 'Seed', seed, 'NormalTransform', 'Inversion'};
else
    % Use ThreeFry for all modern releases as it is faster and better.
    args = {'ThreeFry4x64-20', 'Seed', seed, 'NormalTransform', 'Inversion'};
end
if strcmp( device, 'GPU' )
    rs = parallel.gpu.RandStream(args{:});
else
    rs = RandStream(args{:});
end

%-------------------------------------------------------------------------%
function t = iTimeit( device, f )
% Depending by what device is selected, time f with gputimeit(gpu) or
% timeit(cpu)

% Always disable the "too fast" warnings since we can easily trigger them
% at small sizes. These will be restored to previous state on exit.
warnStruct = [
    warning('off','MATLAB:timeit:HighOverhead')
    warning('off','parallel:gpu:gputimeit:HighOverhead')
    ];
cleaner = onCleanup( @() warning(warnStruct) );

if strcmp( device, 'GPU' )
    t = gputimeit( f );
else
    t = timeit( f );
end

%-------------------------------------------------------------------------%
function iThrowNoGPUError()
error( 'GPUBench:NoGPU', 'No GPU was available for GPUBench to use.' );

%-------------------------------------------------------------------------%
function idx = iMaybeUpdateWaitbar(showWaitbar, idx, numTasks, figh)
idx = idx + 1;
if showWaitbar
    if isvalid(figh)
        waitbar(idx./numTasks, figh);
    else
        % Treat the user closing the progress bar as a request to stop
        error('gpuBench:Cancel', 'User cancelled')
    end
end
