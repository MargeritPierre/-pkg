function mesh = setElementOrder(this,order)
% Change the order of LagrangeElements in the mesh
% Order is a positive or null SCALAR integer  
% Return a copy of the mesh if an output is queried (the original mesh will
% remain unchanged)
% <TODO> how MeshFunctions can react to this to stay valid ?
if nargout==0 ; mesh = this ; else ; mesh = copy(this) ; end

% Get the current element types
    oldElmtTypes = mesh.Elems.Types ;
    nElmtTypes = numel(mesh.Elems.Types) ;

% Depending on the mesh element class
    switch class(oldElmtTypes)
        case 'pkg.geometry.mesh.elements.LagrangeElement'
            % Check if the order shift is possible
                if range([oldElmtTypes.Order])>0 ; error('Order change on mixed-order Lagrange element meshes is not supported') ; end
            % Create the element types of replacement
                newElmtTypes = pkg.geometry.mesh.elements.AbstractElement.empty ;
                for ee = 1:nElmtTypes
                    newElmtTypes(end+1) = pkg.geometry.mesh.elements.LagrangeElement(oldElmtTypes(ee).Geometry,order) ;
                end
        otherwise
            error('The mesh element types are not supported for order change') ;
    end
    
% Set the new element types
    mesh.setElementTypes(newElmtTypes) ;
    

end

