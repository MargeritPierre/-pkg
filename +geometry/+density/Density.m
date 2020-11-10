classdef (Abstract) Density
%DENSITY Abstract class of a density function
    
methods (Abstract)
    h = evalAt(this,p) ;
end

end

