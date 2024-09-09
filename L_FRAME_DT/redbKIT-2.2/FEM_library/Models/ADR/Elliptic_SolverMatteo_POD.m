function [u_FE, MESH, DATA_THERMO] = Elliptic_SolverMatteo_POD(dim, elements, vertices, boundaries, fem, data_file_temp, param, vtk_filename,i_thermo,V_POD,pb_numb)
%ELLIPTIC_SOLVER diffusion-transport-reaction finite element solver
%
%   [U, FE_SPACE, MESH, DATA, ERRORL2, ERRORH1] = ...
%    ELLIPTIC2D_SOLVER(DIM, ELEMENTS, VERTICES, BOUNDARIES, FEM, DATA_FILE, 
%                      PARAM, VTK_FILENAME)
%
%   Inputs:
%     DIM: space dimension, either 2 or 3
%     ELEMENTS, VERTICES, BOUNDARIES: mesh information
%     FEM: string 'P1' or 'P2'
%     DATA_FILE: name of the file defining the problem data and
%          boundary conditions.
%     PARAM: vector of parameters possibly used in the data_file; 
%         if not provided, the PARAM vector is set to the empty vector.
%     VTK_FILENAME: string containing the filename for exporting the
%         solution in the VTK File Format. If not provided or empty, the
%         solution is not exported to vtk.
%
%   Outputs:
%     U: problem solution
%     FE_SPACE: struct containing Finite Element Space information
%     MESH: struct containing mesh information
%     DATA: struct containing problem data

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if isempty(data_file_temp)
    error('Missing data_file_temp')
end

%% Read problem parameters and BCs from data_file
DATA_THERMO   = read_DataFile_temp(data_file_temp, dim, param);
DATA_THERMO.param = param;


%% Set quad_order
if dim == 2
    quad_order       = 4;
elseif dim == 3
    quad_order       = 5;
end

%% Create and fill the MESH data structure
[ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA_THERMO);

%% Create and fill the FE_SPACE data structure
[ FE_SPACE ] = buildFESpace( MESH, fem, 1, quad_order );

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf('-------------------------------------------\n');

%% Preconditioner (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA_THERMO.Preconditioner.type, DATA_THERMO);

%% Assemble matrix and right-hand side
fprintf('\n Assembling ... ');
t_assembly = tic;
[A, F]  =  ADR_Assembler(MESH, DATA_THERMO, FE_SPACE);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

%% Apply boundary conditions
fprintf('\n Apply boundary conditions ');
[A_in, F_in, u_D]   =  ADR_ApplyBC(A, F, FE_SPACE, MESH, DATA_THERMO);

%PROIETTO ALLA GALERKIN
A = V_POD' * A_in * V_POD;
F = V_POD' * F_in;

%% Solve
LinSolver = LinearSolver( DATA_THERMO.LinearSolver );
u_FE = zeros(MESH.numNodes,1);

fprintf('\n Solve Au = f ... ');
Precon.Build( A );
fprintf('\n       **  time to build the preconditioner %3.3f s \n', Precon.GetBuildTime());
LinSolver.SetPreconditioner( Precon );
u = LinSolver.Solve( A, F );
fprintf('\n       ** time to solve the linear system in %3.3f s \n\n', LinSolver.GetSolveTime());

u_FE(MESH.internal_dof) = V_POD * u;
u_FE(MESH.Dirichlet_dof) = u_D;

%% Export to VTK
% if ~isempty(vtk_filename)
%     oldFolder = cd(['/Users/matteotorzoni/Desktop/Corigliano/PORTALE/Dati/portale_calore/Danneggiamento/istanze_',num2str(pb_numb),'/Thermo&Load_',num2str(i_thermo)]);
%     ADR_export_solution(MESH.dim, u_FE(1:MESH.numVertices), MESH.vertices, MESH.elements, vtk_filename);
%     cd(oldFolder);
% end
return
