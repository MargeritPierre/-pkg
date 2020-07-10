classdef MeshFunction < matlab.mixin.Copyable
%MESHFUNCTION Function defined on a mesh
%
    
%% BASE PROPERTIES
properties
    % The parent mesh
    Mesh pkg.geometry.mesh.Mesh
    % Values mesh nodes
    % Size: [nNodes [size(this)]
    Values
end
    
%% CONSTRUCTOR & DESTRUCTOR
methods
    function this = MeshFunction(varargin)
    % Process input arguments
        for arg = 1:2:nargin-1
            this.(varargin{arg}) = varargin{arg+1} ;
        end
    end
    
    function delete(this)
    % Class destructor
    end
end

%% MESH CHANGE HANDLING
properties (Hidden)
    MeshListeners
end
methods
    function setListeners(this)
    % Re-set the mesh listeners
        delete(this.MeshListeners) ;
        this.MeshListeners = event.listener.empty ;
        if ~isempty(this.Mesh)
            this.MeshListeners(end+1) = listener(this.Mesh,'Elems','PostSet',@this.meshChanged) ;
        end
    end
    
    function set.Mesh(this,mesh)
    % Change the mesh
        this.Mesh = mesh ;
        this.setListeners ;
        this.meshChanged ;
    end
    
    function meshChanged(this,src,evt)
    % React to te change of a mesh
%         switch evt.Source{1}.Name
%             case 'Elems' % React to element list change
                % Changes in the node indices
                if this.Mesh.nNodes > size(this.Values,1)
                % New node, add dummy values
                    this.Values(end+1:this.Mesh.nNodes,:) = NaN ;
                elseif this.Mesh.nNodes < size(this.Values,1)
                % Deleted nodes, remove values
                    this.Values(this.Mesh.nNodes+1:end,:) = [] ; % NaN ;
                end
%             otherwise
%         end
    end
end

%% INFOS
methods
    % Size of the function space
    function sz = size(this) ; sz = size(this.Values,2:ndims(this.Values)) ; end
end

%% DATA HANDLING
methods
    function DATA = valuesAtIndices(IND)
    % return data DATA corresponding to the function values at node indices IND
    % See pkg.data.meanDataAtIndices
        DATA = pkg.data.dataAtIndices(this.Values,IND) ;
    end
end

end

