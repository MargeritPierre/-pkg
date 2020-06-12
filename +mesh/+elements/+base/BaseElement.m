classdef (Abstract) BaseElement < pkg.mesh.elements.AbstractElement
%BASEELEMENT Abstract superclass for BASE mesh elements

methods (Sealed) % logical tests sealing
    function [uTypes,ia,ic] = unique(types,varargin)
    % Used to return a unique list of types
        classes = arrayfun(@class,types,'UniformOutput',false) ;
        [~,ia,ic] = unique(classes,varargin{:}) ;
        uTypes = types(ia) ;
    end
    
    function val = eq(a,b)
        classA = arrayfun(@class,a,'UniformOutput',false) ;
        classB = arrayfun(@class,b,'UniformOutput',false) ;
        val = reshape(strcmp(classA,classB),size(a)) ;
    end
end

end

