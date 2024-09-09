function [u, FE_SPACE, MESH, DATA, nodeID, dofID] = ...
    CSMt_PODDEIM_Solver(dim, elements, vertices,...
    boundaries, fem, data_file, param, sol_history, ROM, node_positions,...
    direction_nodeID, instance_time,n_mon_node_out, nodeID, dofID, i3, index_param)
%CSMT_PODDEIM_SOLVER Dynamic Structural POD-DEIM Solver

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 6
    error('Missing input arguments. Please type help CSM_Solver')
end

if isempty(data_file)
    error('Missing data_file')
end

if nargin < 7
    param = [];
end

if nargin < 8
    vtk_filename = [];
end

if nargin < 9 || isempty(sol_history)
    sol_history = false;
end

% counter_vtk_out       = n_vtk_out;
% counter_done_out      = 0;
counter_mon_node_out  = n_mon_node_out;
counter_node_done_out = 0;

%% Read problem parameters and BCs from data_file
DATA   = CSM_read_DataFile(data_file, dim, param);
if nargin < 7
    DATA.param = [];
else
    DATA.param = param;
end
t      = [];

%% Set quad_order
if dim == 2
    quad_order       = 4;
elseif dim == 3
    quad_order       = 5;
end

%% Create and fill the MESH data structure
[ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA, 'CSM' );

%% Gather Time Setting
t0        = DATA.time.t0;
dt        = DATA.time.dt;
if instance_time
    tf        = DATA.time.tf_long;
else
    tf        = DATA.time.tf_short;
end
t         = DATA.time.t0;
k_t       = 0;

TimeAdvance = GeneralizedAlpha_TimeAdvance( DATA.time.beta, DATA.time.gamma, DATA.time.alpha_m, DATA.time.alpha_f, dt );

%% Create and fill the FE_SPACE data structure
[ FE_SPACE ] = buildFESpace( MESH, fem, dim, quad_order );

%% Determine the dofs to be saved
if ~isempty(n_mon_node_out) && isempty(dofID)
    % estimate the number of required steps
    n_steps=floor(tf/dt/n_mon_node_out);
    % find the ID of the monitored nodes
    nodeID=ones(1,size(node_positions,2));
    for i1=1:size(node_positions,2)
        translated_vertices = MESH.vertices - node_positions(:,i1);
        norm_translation=100.*ones(1,size(translated_vertices,2));
        for i2 = 1:length(norm_translation)
            if MESH.dim == 2
                norm_translation(i2) = (translated_vertices(1,i2))^2+(translated_vertices(2,i2))^2;
            elseif MESH.dim == 3
                norm_translation(i2) = (translated_vertices(1,i2))^2+(translated_vertices(2,i2))^2+(translated_vertices(3,i2))^2;
            end
        end
        [~,nodeID(i1)]=min(norm_translation);
    end
    dofID=zeros(1,size(node_positions,2));
    for i1=1:length(nodeID)
        dofID(i1) = nodeID(i1)+(direction_nodeID(i1))*MESH.numNodes;
    end
    monitored_gdls.time_steps=zeros(n_steps,1);
    monitored_gdls.displ=zeros(n_steps,length(nodeID));
    
elseif isempty(n_mon_node_out)
    dofID=[];
    nodeID=[];
    
else
    % estimate the number of required steps
    n_steps=floor(tf/dt/n_mon_node_out);
    monitored_gdls.time_steps=zeros(n_steps,1);
    monitored_gdls.displ=zeros(n_steps,length(nodeID));
    
end


u0  = [];
du0 = [];
for k = 1 : FE_SPACE.numComponents
    switch dim
        case 2
            u0  = [u0; DATA.u0{k}(  MESH.nodes(1,:), MESH.nodes(2,:), t0, param )'];
            du0 = [du0; DATA.du0{k}( MESH.nodes(1,:), MESH.nodes(2,:), t0, param )'];
            
        case 3
            u0  = [u0; DATA.u0{k}(  MESH.nodes(1,:), MESH.nodes(2,:), MESH.nodes(3,:), t0, param )'];
            du0 = [du0; DATA.du0{k}( MESH.nodes(1,:), MESH.nodes(2,:), MESH.nodes(3,:), t0, param )'];
    end
end

u = u0;
% if ~isempty(vtk_filename)
%     CSM_export_solution(MESH.dim, u0, MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, 0);
% end

Coef_Mass = TimeAdvance.MassCoefficient( );

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices         = %d \n', MESH.numVertices);
fprintf(' * Number of Nodes            = %d \n', MESH.numNodes);
fprintf(' * Number of Elements         = %d \n', MESH.numElem);
fprintf(' * Number of Reduced Elements = %d \n', ROM.Red_Mesh.numElem);
fprintf(' * %% Reduced Elements         = %2.2f \n', ROM.Red_Mesh.numElem / MESH.numElem * 100);
fprintf(' * Number of Dofs             = %d \n', size(ROM.V,2));
fprintf(' * Number of timesteps        = %d \n', (tf-t0)/dt);
fprintf('-------------------------------------------\n');

%% Preconditioner (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Preconditioner.type, DATA);

SolidModel = CSM_Assembler( ROM.Red_Mesh, DATA, FE_SPACE );

%% Assemble mass matrix
fprintf('\n Assembling mass matrix... ');
t_assembly = tic;
M    =  DATA.Density * ROM.M;
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

LinSolver = LinearSolver( DATA.LinearSolver );

%% Initial Acceleration
fprintf('\n -- Assembling external Forces at t0... ');
t_assembly = tic;
F_ext_0_FE = SolidModel.compute_volumetric_forces(t0);
F_ext_0_FE = F_ext_0_FE( MESH.internal_dof );
F_ext_0    = ROM.LeftProjection_ext * F_ext_0_FE( ROM.IDEIM_ext );
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

fprintf('\n -- Assembling internal Forces at t0... ');
t_assembly = tic;
F_in_0_FE  = SolidModel.compute_internal_forces( u0, t );
F_in_0_FE  = F_in_0_FE( MESH.internal_dof );
F_in_0     = ROM.LeftProjection_int * F_in_0_FE( ROM.IDEIM_in );
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly)

d2u0_N = M \ (F_ext_0 - F_in_0);

TimeAdvance.Initialize( ROM.V'*u0(MESH.internal_dof), ROM.V'*du0(MESH.internal_dof), d2u0_N);%ROM.V'*d2u0(MESH.internal_dof) );

U_kN = ROM.V' * u(MESH.internal_dof,1);
U_n  = u0;

%% Time Loop
while ( t < tf )
    
    iter_time = tic;
    
    t       = t   + dt;
    k_t     = k_t + 1;
    
    fprintf('\n=========================================================================')
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n', t0, t, tf );

    % Newton Method
    tol        = DATA.NonLinearSolver.tol;
    resRelNorm = tol + 1;
    incrNorm   = tol + 1;
    maxIter    = DATA.NonLinearSolver.maxit;
    k          = 1;
    
    [~, ~, u_D]   =  CSM_ApplyEssentialBC([], [], MESH, DATA, t);
    dU             = zeros(MESH.numNodes*MESH.dim,1);
    U_k            = u(:,end);
    U_k(MESH.Dirichlet_dof) = u_D;
    
    Csi = TimeAdvance.RhsContribute( );
    
    % Assemble matrix and right-hand side
    fprintf('\n -- Assembling external Forces... ');
    t_assembly = tic;
    F_ext_FE      = SolidModel.compute_external_forces( (1 - TimeAdvance.M_alpha_f) * t + TimeAdvance.M_alpha_f * (t-dt) );
    F_ext_FE      = F_ext_FE( MESH.internal_dof );
    F_ext         = ROM.LeftProjection_ext * F_ext_FE( ROM.IDEIM_ext );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    fprintf('\n -- Assembling internal Forces ... ');
    t_assembly = tic;
    F_in_FE      = SolidModel.compute_internal_forces( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n, t );
    F_in_FE      = F_in_FE(MESH.internal_dof);
    F_in         = ROM.LeftProjection_int * F_in_FE( ROM.IDEIM_in );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    Residual  = Coef_Mass * M * U_kN + F_in - F_ext - M * Csi;
        
    fprintf('\n -- Assembling Jacobian matrix... ');
    t_assembly = tic;
    dF_in_FE   = SolidModel.compute_jacobian( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n, t );
    dF_in_FE   = CSM_ApplyEssentialBC(dF_in_FE, [], MESH, DATA, t, 1);
    dF_in      = ROM.LeftProjection_int * ( dF_in_FE(ROM.IDEIM_in, : ) * ROM.V );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);        
            
    Jacobian  = Coef_Mass * M + (1 - TimeAdvance.M_alpha_f) * dF_in;

    res0Norm = norm(Residual);
        
    fprintf('\n============ Start Newton Iterations ============\n\n');
    while (k <= maxIter && (incrNorm > tol || resRelNorm > tol) )
        
        % Solve
        fprintf('\n   -- Solve J x = -R ... ');
        Precon.Build( Jacobian );
        fprintf('\n        time to build the preconditioner %3.3f s \n', Precon.GetBuildTime());
        LinSolver.SetPreconditioner( Precon );
        dU_N                  = LinSolver.Solve( Jacobian, -Residual );
        dU(MESH.internal_dof) = ROM.V * dU_N;
        fprintf('\n        time to solve the linear system in %3.3f s \n', LinSolver.GetSolveTime());
        
        % update solution
        U_k        = U_k + dU;
        U_kN       = U_kN + dU_N;
        incrNorm   = norm(dU)/norm(U_k);
        
        % Assemble matrix and right-hand side
        fprintf('\n   -- Assembling internal forces... ');
        t_assembly = tic;
        F_in_FE      = SolidModel.compute_internal_forces( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n, t );
        F_in_FE      = F_in_FE(MESH.internal_dof);
        F_in         = ROM.LeftProjection_int * F_in_FE( ROM.IDEIM_in );
        dF_in_FE     = SolidModel.compute_jacobian( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n, t );
        dF_in_FE     = CSM_ApplyEssentialBC(dF_in_FE, [], MESH, DATA, t, 1);
        dF_in        = ROM.LeftProjection_int * ( dF_in_FE(ROM.IDEIM_in, : ) * ROM.V );
        t_assembly = toc(t_assembly);
        fprintf('done in %3.3f s\n', t_assembly);
        
        Residual  = Coef_Mass * M * U_kN + F_in - F_ext - M * Csi;
    
        Jacobian  = Coef_Mass * M + (1 - TimeAdvance.M_alpha_f) * dF_in;
        
        resRelNorm = norm(Residual) / res0Norm;
        
        fprintf('\n **** Iteration  k = %d:  norm(dU)/norm(Uk) = %1.2e, Residual Rel Norm = %1.2e \n\n',k,full(incrNorm), full(norm(resRelNorm)));
        k = k + 1;
        
    end
    
    if sol_history
        u = [u U_k];
    else
        u = U_k;
    end
    
    if ~isempty(counter_mon_node_out)
        if counter_mon_node_out==0
            counter_node_done_out=counter_node_done_out+1;
            monitored_gdls.time_steps(counter_node_done_out) = t;
            monitored_gdls.displ(counter_node_done_out,:)= U_n(dofID);
            counter_mon_node_out=n_mon_node_out;
        else
            counter_mon_node_out=counter_mon_node_out-1;
        end
    end
    
    U_n = U_k;

    %% Export to VTK
%     if ~isempty(vtk_filename)
%         CSM_export_solution(MESH.dim, U_k, MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, k_t);
%     end
    
    TimeAdvance.Update( U_kN );
    
    iter_time = toc(iter_time);
    fprintf('\n-------------- Iteration time: %3.2f s -----------------',iter_time);
    
end

if ~isempty(dofID)
    number_entries=find(monitored_gdls.time_steps==0)-1;
    che = monitored_gdls.time_steps - DATA.time.t_in_record;
    [~,guevara] = min(abs(che));
    
    monitored_gdls.time_steps=monitored_gdls.time_steps(guevara:number_entries);
    monitored_gdls.displ=monitored_gdls.displ(guevara:number_entries,:);
   
    if ~isempty(index_param)
        [~, ~, ~] = mkdir(['C:\Users\Rosafalco\Documents\Polimi\3_Dottorato\Dati\Non_linear_fitting\Damage_Identification\mensola\damaged_',DATA.pb_ID_f,'\damaged_', DATA.pb_ID], 'graph_check');
        save(['C:\Users\Rosafalco\Documents\Polimi\3_Dottorato\Dati\Non_linear_fitting\Damage_Identification\mensola\damaged_', DATA.pb_ID_f, '\damaged_', DATA.pb_ID, '\graph_check\monitored_gdls_U_POD_DEIM', num2str(index_param),'.mat']','-struct', 'monitored_gdls');
        save(['C:\Users\Rosafalco\Documents\Polimi\3_Dottorato\Dati\Non_linear_fitting\Damage_Identification\mensola\damaged_', DATA.pb_ID_f, '\damaged_', DATA.pb_ID, '\damage_level_1\monitored_gdls_U', num2str(i3),'.mat']','-struct', 'monitored_gdls');
    end
    save(['C:\Users\Rosafalco\Documents\Polimi\3_Dottorato\Dati\Non_linear_fitting\Damage_Identification\mensola\damaged_', DATA.pb_ID_f, '\damaged_', DATA.pb_ID, '\damage_level_1\monitored_gdls_U', num2str(i3),'.mat']','-struct', 'monitored_gdls');
end

fprintf('\n************************************************************************* \n');

return
