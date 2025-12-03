function plotFLOPS( data, color, linestyle, linewidth )
%plotFLOPS  plot some gpuBench data as floating-point operations per second
%
%   gpubench.plotFLOPS(data,color,linewidth)

%   Copyright 2011-2024 The MathWorks, Inc.

semilogx(data.Sizes, 1e-12 * double(data.NumOps) ./ data.Times, ...
    'Color', color, ...
    'Marker', '.', ...
    'MarkerSize', 10+6*linewidth, ...
    'linestyle', linestyle, ...
    'linewidth', linewidth)
xlabel( 'Number of elements' )
ylabel( sprintf('TFLOPS\n(higher is better)') )
grid( 'on' )
hold( 'on' )

% Fix the position to avoid over-size images on high-DPI displays
pos = get(gcf, 'Position');
pos(:,3:4) = [650 400];
set(gcf, 'Position',pos)