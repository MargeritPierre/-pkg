classdef ElementMesh
%ELEMENTMESH A try to implement a mesh with elements
    
properties
    Elements = pkg.mesh.elements.ElementTable ;
    Faces = pkg.mesh.elements.ElementTable ;
    Edges = pkg.mesh.elements.ElementTable ;
end
    
methods
    function obj = ElementMesh(varargin)
    % Class constructor
    end
end
end

