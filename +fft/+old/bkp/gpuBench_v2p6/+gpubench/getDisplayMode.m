function [showProgress,canShowGraphics] = getDisplayMode()
%getDisplayMode  determine whether and how to indicate progress.

%   Copyright 2025 The MathWorks, Inc.

% For now hard-code that we always show progress. May become an option in
% future.
showProgress = true;

% Check if graphical display is enabled.
if verLessThan('MATLAB', '9.9') %#ok<VERLESSMATLAB> before R2020b
    % Old versions do not provide an easy way to check display availability, so
    % just go for the most common value.
    canShowGraphics = true;
elseif isMATLABReleaseOlderThan('R2024a')
    canShowGraphics = iGetDisplayCapabilityOld();
else
    canShowGraphics = iGetDisplayCapabilityLatest();
end

end

function tf = iGetDisplayCapabilityLatest()
import matlab.internal.capability.*;
tf = Capability.isSupported(Capability.ModalDialogs);
end

function tf = iGetDisplayCapabilityOld()
import matlab.internal.lang.capability.*;
tf = Capability.isSupported(Capability.ModalDialogs);
end
