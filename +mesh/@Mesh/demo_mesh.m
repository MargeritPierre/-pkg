%% TEST MESH PACKAGE
clc,clearvars
import pkg.mesh.*
import pkg.mesh.elements.*


% mesh = Mesh('Elems',[1 2 3 NaN ; 3 1 4 5]) ;
% 
% mesh.Elems = [mesh.Elems [3 4 6]] ; % element table concatenation
% mesh.X = rand(mesh.nNodes,3) ; % set mesh nodal coordinates

%
nNodes = 1000000 ;
nElems = 100000 ;
elemType = LagrangeElement('quad',1) ;
nodeIndices = randi([1 nNodes],[nElems elemType.nNodes]) ;

elemTable = ElementTable('Types',elemType,'Indices',nodeIndices) ;
mesh = Mesh('Elems',elemTable) ;

tic ;
profile on

edges = mesh.Elems.allEdges ;

profile viewer
toc

