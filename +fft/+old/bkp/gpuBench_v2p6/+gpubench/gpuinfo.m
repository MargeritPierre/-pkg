function info = gpuinfo()
%GPUINFO  get basic details about the GPU being used
%
%   GPUINFO() gets hold of some basic info about the GPU device being used.
%
%   See also: CPUINFO

%   Copyright 2011-2025 The MathWorks, Inc.

if parallel.gpu.GPUDevice.isAvailable()
    gpu = gpuDevice();
    info = struct( ...
        'Name', gpu.Name, ...
        'Clock', sprintf( '%u MHz', gpu.ClockRateKHz/1e3 ), ...
        'NumProcessors', gpu.MultiprocessorCount, ...
        'ComputeCapability', gpu.ComputeCapability, ...
        'TotalMemory', sprintf( '%1.2f GB', gpu.TotalMemory/2^30 ), ...
        'CUDAVersion', gpu.DriverVersion, ...
        'DriverVersion', iGetGraphicsDriverVersion(gpu) );
else
    % No GPU, so create an empty structure
    info = struct( ...
        'Name', '<no GPU available>', ...
        'Clock', '', ...
        'NumProcessors', 0, ...
        'ComputeCapability', '', ...
        'TotalMemory', 0, ...
        'CUDAVersion', '', ...
        'DriverVersion', '' );
end

end

function v = iGetGraphicsDriverVersion(gpu)
% Try to determine the grpahics driver version. This is a property of
% gpuDevice for releases R2023a and newer, but for older releases we try
% the old way and fallback to "unknown" if that doesn't work.
if isprop(gpu,"GraphicsDriverVersion")
    v = gpu.GraphicsDriverVersion;
else
    % Make sure we don't error if the driver string could not be found.
    try
        v = parallel.internal.gpu.CUDADriverVersion();
    catch
        v = "unknown";
    end
end
end