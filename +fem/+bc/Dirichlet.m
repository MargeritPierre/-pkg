classdef Dirichlet < pkg.fem.bc.AbstractBC
%DIRICHLET boundary condition
    
properties
    Near % function handle @(x) or coordinates, where to apply the BC
    Values % function handle @(x) or vector, BC value
end

methods
    function this = Dirichlet(near,val)
    % Class constructor
        if nargin<2 ; error('Two arguments must be provided') ; end
        this.Near = near ;
        this.Values = val(:)' ;
    end

    function [v0,fixDOF] = apply(this,mesh)
    % Everything needed to apply the BC to a mesh
        fixDOF = mesh.near(this.Near) & mesh.boundaryNodes ;
        if isa(this.Values,'function_handle') 
            vBC = this.Values(mesh.Nodes(fixDOF,:)) ;
        else
            vBC = repmat(this.Values,sum(fixDOF),1) ;
        end
        if nargout<2
            v0 = NaN(mesh.nNodes,size(vBC,2)) ;
        else
            v0 = zeros(mesh.nNodes,size(vBC,2)) ;
        end
        v0(fixDOF,:) = vBC ;
    end
end

end

