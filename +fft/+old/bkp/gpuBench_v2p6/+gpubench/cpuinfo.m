function info = cpuinfo()
%CPUINFO  read CPU configuration
%
%   info = CPUINFO() returns a structure containing various bits of
%   information about the CPU and operating system as provided by /proc/cpu
%   (Unix), sysctl (Mac) or WMIC (Windows). This information includes:
%     * CPU name
%     * CPU clock speed
%     * CPU Cache size (L2)
%     * Number of physical CPU cores
%     * Operating system name & version
%
%   See also: COMPUTER, ISUNIX, ISMAC

%   Author: Ben Tordoff
%   Copyright 2011-2025 The MathWorks, Inc.

if isunix
    if ismac
        info = cpuInfoMac();
    else
        info = cpuInfoUnix();
    end
else
    info = cpuInfoWindows();
end


%-------------------------------------------------------------------------%
function info = cpuInfoWindows()
sysInfo = callWMI( 'Processor' ); 
osInfo = callWMI( 'OperatingSystem' );

info = struct( ...
    'Name', strtrim(char(sysInfo(1).Name)), ...
    'Clock', [num2str(sysInfo(1).MaxClockSpeed),' MHz'], ...
    'Cache', [num2str(sysInfo(1).L2CacheSize),' KB'], ...
    'NumProcessors', sum( [sysInfo.NumberOfCores] ), ...
    'OSType', 'Windows', ...
    'OSVersion', strtrim(char(osInfo.Caption)) );

%-------------------------------------------------------------------------%
function info = callWMI( alias )
% Call the MS-DOS WMI (Windows Management) command

NET.addAssembly('System.Management');
objects = System.Management.ManagementObjectSearcher("SELECT * FROM Win32_" + alias).Get;
numResults = objects.Count;

% Get the returned object and convert its contents into a struct.
objEnum = objects.GetEnumerator;

objEnum.MoveNext;
object = objEnum.Current;
props = object.Properties;
N = props.Count;
fields = cell(N,1);
values = cell(N,numResults);
for jj=1:numResults
    ii = 0;
    e = props.GetEnumerator;
    while e.MoveNext
        ii = ii + 1;
        tmp = e.Current;
        fields{ii} = char(tmp.Name);
        values{ii,jj} = tmp.Value;
    end
    objEnum.MoveNext;
end

% Convert to a structure
info = cell2struct( values, fields );

%-------------------------------------------------------------------------%
function info = cpuInfoMac()
machdep = callSysCtl( 'machdep.cpu' );
hw = callSysCtl( 'hw' );
if strcmp(computer,'MACA64')
    maxFreq = 'N/A';
    cache = 'N/A';
else
    maxFreq = [num2str(str2double(hw.cpufrequency_max)/1e6),' MHz'];
    cache = [machdep.cache.size,' KB'];
end

info = struct( ...
    'Name', machdep.brand_string, ...
    'Clock', maxFreq, ...
    'Cache', cache, ...
    'NumProcessors', str2double( machdep.core_count ), ...
    'OSType', 'Mac OS/X', ...
    'OSVersion', getOSXVersion() );

%-------------------------------------------------------------------------%
function info = callSysCtl( namespace )
[~, infostr] = system( sprintf( 'sysctl -a %s', namespace ) );
% Remove the prefix
infostr = strrep( infostr, [namespace,'.'], '' );
% Now break into a structure
infostr = textscan( infostr, '%s', 'delimiter', '\n' );
infostr = infostr{1};
info = struct();
for ii=1:numel( infostr )
    colonIdx = find( infostr{ii}==':', 1, 'first' );
    if isempty( colonIdx ) || colonIdx==1 || colonIdx==length(infostr{ii})
        continue
    end
    prefix = infostr{ii}(1:colonIdx-1);
    suffix = '';
    value = strtrim(infostr{ii}(colonIdx+1:end));
    while ismember( '.', prefix )
        dotIndex = find( prefix=='.', 1, 'last' );
        suffix = prefix(dotIndex+1:end);
        prefix = prefix(1:dotIndex-1);
    end
    if ~isempty(suffix)
        info.(prefix).(suffix) = value;
    else
        info.(prefix) = value;
    end    
end

%-------------------------------------------------------------------------%
function vernum = getOSXVersion()
% Extract the OS version number from the system software version output.
[~, ver] = system('sw_vers');
vernum = regexp(ver, 'ProductVersion:\s([1234567890.]*)', 'tokens', 'once');
vernum = strtrim(vernum{1});

%-------------------------------------------------------------------------%
function info = cpuInfoUnix()
txt = readCPUInfo();
cpuinfo = parseCPUInfoText( txt );

txt = readOSInfo();
osinfo = parseOSInfoText( txt );

% Merge the structures
info = cell2struct( [struct2cell( cpuinfo );struct2cell( osinfo )], ...
    [fieldnames( cpuinfo );fieldnames( osinfo )] );

%-------------------------------------------------------------------------%
function info = parseCPUInfoText( txt )
% Now parse the fields
lookup = {
    'model name', 'Name'
    'cpu Mhz', 'Clock'
    'cpu cores', 'NumProcessors'
    'cache size', 'Cache'
    };
info = struct( ...
    'Name', {''}, ...
    'Clock', {''}, ...
    'Cache', {''}, ...
    'NumProcessors', {[]} );
for ii=1:numel( txt )
    if isempty( txt{ii} )
        continue;
    end
    % Look for the colon that separates the property name from the value
    colon = find( txt{ii}==':', 1, 'first' );
    if isempty( colon ) || colon==1 || colon==length( txt{ii} )
        continue;
    end
    fieldName = strtrim( txt{ii}(1:colon-1) );
    fieldValue = strtrim( txt{ii}(colon+1:end) );
    if isempty( fieldName ) || isempty( fieldValue )
        continue;
    end
    
    % Is it one of the fields we're interested in?
    idx = find( strcmpi( lookup(:,1), fieldName ) );
    if ~isempty( idx )
        newName = lookup{idx,2};
        info.(newName) = fieldValue;
    end
end

% Convert clock speed
info.Clock = [info.Clock, ' MHz'];

% Convert num cores
info.NumProcessors = str2double( info.NumProcessors );

%-------------------------------------------------------------------------%
function info = parseOSInfoText( txt )
info = struct( ...
    'OSType', 'Linux', ...
    'OSVersion', '' );
% find the string "linux version" then look for the bit in brackets
[~,b] = regexp( txt, '[^\(]*\(([^\)]*)\).*', 'match', 'tokens', 'once' );
info.OSVersion = b{1}{1};

%-------------------------------------------------------------------------%
function txt = readCPUInfo()

fid = fopen( '/proc/cpuinfo', 'rt' );
if fid<0
    error( 'cpuinfo:BadPROCCPUInfo', 'Could not open /proc/cpuinfo for reading' );
end
cleanup = onCleanup( @() fclose( fid ) );

txt = textscan( fid, '%s', 'Delimiter', '\n' );
txt = txt{1};

%-------------------------------------------------------------------------%
function txt = readOSInfo()

fid = fopen( '/proc/version', 'rt' );
if fid<0
    error( 'cpuinfo:BadProcVersion', 'Could not open /proc/version for reading' );
end
cleanup = onCleanup( @() fclose( fid ) );

txt = textscan( fid, '%s', 'Delimiter', '\n' );
txt = txt{1};
