function setPlot(gammaMax,plotType,logScale)
% SET A WAVE DISPERSION PLOT GIVEN SOME INPUT OPTIONS

    % Get lines
    lines = findobj(gca,'type','line')' ;
    lines = [lines findobj(gca,'type','scatter')'] ;
    lines = [lines findobj(gca,'type','patch')'] ;

    % CullgGamma>gammaMax
    for ll = lines 
        ki = ll.YData + 1i*ll.ZData ;
        valid = ones(size(ki)) ; 
        valid(abs(imag(ki)./real(ki))>gammaMax) = NaN ; 
        ll.YData = ll.YData.*valid ;
        ll.ZData = ll.ZData.*valid ;
    end

    % Wavenumber, phase velocity, etc
    xlabel 'Frequency (Hz)'
    switch plotType
        case 'c'
            ylabel 'Phase Velocity (m/s)'
            for ll = lines 
                omega = 2*pi*ll.XData ;
                ki = ll.YData + 1i*ll.ZData ;
                ll.YData = 1e-3*real(omega./ki) ; 
                ll.ZData = 1e-3*imag(omega./ki) ; 
            end
        case 'k'
            ylabel 'Wavenumber (rad/mm)'
        case 'g'
            ylabel 'Wavenumber (rad/mm)'
            zlabel 'Spatial decay'
            for ll = lines 
                ki = ll.YData + 1i*ll.ZData ;
                ll.ZData = imag(ki)./real(ki) ; 
            end
        case 'l'
            ylabel 'Wavelength (mm)'
            for ll = lines 
                ki = ll.YData + 1i*ll.ZData ;
                ll.YData = real(2*pi./ki) ; 
                ll.ZData = imag(2*pi./ki) ; 
            end
    end

    % Scale change
    if logScale
        for ll = lines 
            ll.YData = abs(ll.YData) ; 
            ll.ZData = abs(ll.ZData) ; 
        end
        set(gca,'xscale','log','yscale','log','zscale','log')
    end

end
