classdef MeshSpace < pkg.geometry.mesh.elements.TableObject
%MESHSPACE Support for a mesh function defining the interpolation of a
%function in the mesh
    
properties (SetAccess = immutable)
    Mesh pkg.geometry.mesh.Mesh
    ElemTypes pkg.geometry.mesh.elements.AbstractElement
    MeshListeners event.listener
end

methods
    function this = MeshSpace(mesh,elemTypes)
    % Class constructor
        if nargin<1 ; error('A mesh MUST be provided !') ; end
    % Set the mesh
        this.Mesh = mesh ;
    % Set element types
        if nargin<2 ; this.ElemTypes = this.Mesh.Elems.Types ;
        else ; this.ElemTypes = elemTypes ;
        end
    % Compute the mesh space data
        this.reset ;
    % Listen to mesh changes
        %this.MeshListeners(end+1) = addlistener(this.Mesh,'Nodes','PostSet',@this.reset) ;
        this.MeshListeners(end+1) = addlistener(this.Mesh,'Elems','PostSet',@this.reset) ;
    end
    
    function reset(this,~,~)
    % Reset the mesh space data
        intMesh = this.Mesh.setElementTypes(this.ElemTypes) ;
        this.Elems = intMesh.Elems ;
    end
    
    function delete(this)
    % Class destructor
    end
end

end

