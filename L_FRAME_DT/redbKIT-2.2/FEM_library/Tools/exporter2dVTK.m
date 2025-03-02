function exporter2dVTK(data)
%EXPORTER2DVTK writes a VTK file for unstructured mesh
%
%  Description:
%
%    VTK exporter for 2D TRIANGULAR UNSTRUCTURED FEM simulations.
%    Use Paraview to plot results
%    (adapted from the original work of John Burkardt)
%
%  Author:
%
%    Matteo Astorino (ASCII version)
%    Federico Negri  (BINARY version)
%
%  Parameters:
%
%    coor(nodeNum,2) are the node coordinates.
%
%    elementConnectivity(elementNum,nodePerElement) describe the nodes that
%    form each element.
%
%    U(nodeNum,2), P(1,nodeNum), U contains velocity components and P the
%    pressure at each node.
%
%    OUTPUT_FILENAME, the name of the output file.
%    By convention, this file should have the extension ".vtk".
%
%    TITLE is the title for the data.
%
timewrite = tic;
%
%  Get the size of your mesh
%
[ nodeNum, dim ] = size ( data.vertices );
[ elementNum,  nodePerElement ] = size ( data.elements );
%
%  Open the output file.
%
if data.iteration == -1
    output_filename = sprintf('%s.vtk', data.outputFile);
else
    output_filename = sprintf('%s%04d.vtk', data.outputFile,data.iteration);
end

if ( isempty ( output_filename ) )
    output_filename = 'solution.vtk';
end

outputSolution = fopen ( output_filename, 'w' );
%
%  Transpose or otherwise modify the data.
%
coorz = zeros ( 3, nodeNum );
if size(data.vertices,2) == 3
    coorz(1:3,1:nodeNum) = (data.vertices(1:nodeNum,1:3))';
else
    %L: nel caso di un'analisi bidimensionale si inizializza a 0 la riga
    %riferita alla coordinata z
    coorz(1:2,1:nodeNum) = (data.vertices(1:nodeNum,1:2))';
end


elementNode = zeros ( nodePerElement, elementNum );
%L: la numerazione in PARAVIEW inizia da 0: per questo � necessario
%rieffettuare la numerazione dei nodi
elementNode(1:nodePerElement,1:elementNum) = ...
    data.elements(1:elementNum,1:nodePerElement)' - 1;


%
%  Write the data.
%
vtk_headerfile ( outputSolution, data.title, nodeNum, elementNum, ...
    nodePerElement, coorz, elementNode, output_filename);
for i=1:size(data.variableName,2)
    vtk_write ( outputSolution, nodeNum, data.variableName{i}, data.variableType{i}, data.variableData{i}, output_filename);
end
fclose ( outputSolution );

timewrite = toc(timewrite);
fprintf ( 1, '\n' );
fprintf ( 1, '  The data was written to "%s" in %1.3f seconds\n', output_filename, timewrite);

return
end

function vtk_headerfile ( outputSolution, title, nodeNum, elementNum, ...
    nodePerElement, coorz, elementNode, filename)
%*****************************************************************************
%
%% VTK_HEADERFILE write the header of the VTK file.
%
%  Author:
%
%    Matteo Astorino
%
%  Parameters:
%
%    Input, integer outputSolution, the output unit.
%
%    Input, string TITLE, a title for the data.
%
%    Input, integer nodeNum, the number of nodes.
%
%    Input, integer elementNum, the number of elements.
%
%    Input, integer nodePerElement, the order of the elements.
%
%    Input, real coorZ(3,nodeNum), the node coordinates.
%
%    Input, integer elementNode(nodePerElement,elementNum), the
%    nodes that make up each element.  Node indices are zero-based.
%
%
fprintf ( outputSolution, '# vtk DataFile Version 3.0\n' );
fprintf ( outputSolution, '%s\n', title );
fprintf ( outputSolution, 'BINARY\n' );
fprintf ( outputSolution, '\n' );
fprintf ( outputSolution, 'DATASET UNSTRUCTURED_GRID\n' );

fprintf ( outputSolution, 'POINTS %d float\n', nodeNum );
fwrite(outputSolution, coorz(1:3,:),'float','b');

%
%  Note that CELL_SIZE uses nodePerElement+1 because the order of each element
%  is included as a data item.
%
cell_size = elementNum * ( nodePerElement + 1 );

fprintf ( outputSolution, '\n' );
fprintf ( outputSolution, 'CELLS  %d  %d\n', elementNum, cell_size );

fwrite(outputSolution, [nodePerElement*ones(1,elementNum); elementNode],'int','b');

fprintf ( outputSolution, '\n' );
fprintf ( outputSolution, 'CELL_TYPES %d\n', elementNum );

fwrite(outputSolution, [5*ones(elementNum,1)]','int','b');


fprintf ( outputSolution, '\n' );
fprintf ( outputSolution, 'POINT_DATA %d\n', nodeNum );
return
end

function vtk_write ( outputSolution, nodeNum, variableName, variableType, variableData, filename)
%*****************************************************************************80
%
%% VTK_WRITE vector or scalar data to a VTK file.
%
%  Author:
%
%    Matteo Astorino
%
%  Parameters:
%
%    Input, integer outputSolution, the output unit.
%
%    Input, integer nodeNum, the number of nodes.
%
%    Input, real P(1,nodeNum), the pressure at each node.
%
%    Input, real UVW(3,nodeNum), the velocity at each node.
if(variableType=='SCALARS')
    fprintf ( outputSolution, 'SCALARS %s float\n', variableName );
    fprintf ( outputSolution, 'LOOKUP_TABLE default\n' );
    
    fwrite(outputSolution, variableData(1,:),'float','b');
    fprintf ( outputSolution, '\n' );
    
else
    
    fprintf ( outputSolution, 'VECTORS %s float\n',variableName );
    fwrite(outputSolution, variableData(1:3,:),'float','b');
    fprintf ( outputSolution, '\n' );
    
    return
end
end
