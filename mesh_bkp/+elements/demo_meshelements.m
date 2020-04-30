
%% Import the package
import pkg.mesh.elements.*

% Create the element (choose any class you want)
clc
%elmt = CustomElement ;
elmt = LagrangeElement('hex',2)

% Display the reference element
clf
plot(elmt,'ReferenceElement') ; % equivalent to plot(elmt)
axis equal

%% PLOT THE SHAPE FUNCTIONS
clf
plot(elmt,'ShapeFunctions') ;


%% ELEMENT TABLES
clearvars,clc
import pkg.mesh.elements.*

ElementTypes = LagrangeElement('quad',1) ;
nodeIndices = [1 2 3 4 ; 4 5 6 7] ;

T = ElementTable('NodeIndices',nodeIndices,'Types',ElementTypes) ;









