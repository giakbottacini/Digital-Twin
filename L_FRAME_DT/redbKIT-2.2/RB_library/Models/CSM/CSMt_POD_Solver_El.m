function [u, FE_SPACE, MESH, DATA, nodeID, dofID] = ...
                CSMt_POD_Solver_El(FOM_affine, MESH, FE_SPACE, dim, fem, ...
        data_file, param, node_positions, direction_nodeID, vtk_filename, ...
        Training_Options, instance_time, n_vtk_out, n_mon_node_out, ...
        nodeID, dofID, i3, index_param, V_POD)

%CSMT_POD_SOLVER Elasto-Dynamic Structural POD-Galerkin Solver

%   V_POD is a POD basis ( i.e. matrix of size #InternalDoFs x #PODmodes )

if isempty(data_file)
    error('Missing data_file')
end

if isempty( Training_Options )
    export_h5 = false;
else
    export_h5 = true;
end

% Inizializzazione dei contatori
counter_vtk_out       = n_vtk_out;
counter_done_out      = 0;
counter_mon_node_out  = n_mon_node_out;
counter_node_done_out = 0;

%% Read problem parameters and BCs from data_file
DATA       = CSM_read_DataFile(data_file, dim, param);
DATA.param = param;


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

%% Variable inizialization
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
if ~isempty(vtk_filename)
    CSM_export_solution(MESH.dim, u0, MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, 0);
end

Coef_Mass = TimeAdvance.MassCoefficient( );

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf(' * Number of Dofs      = %d \n',size(V_POD, 2));
fprintf(' * Number of timesteps =  %d\n', (tf-t0)/dt);
fprintf('-------------------------------------------\n');

%% Preconditioner (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Preconditioner.type, DATA);

SolidModel    = CSM_Assembler( MESH, DATA, FE_SPACE );

%% Assemble mass matrix
% fprintf('\n Assembling mass matrix... ');
% t_assembly = tic;
M_h   = SolidModel.compute_mass();
M_h   = M_h * DATA.Density;
M     = V_POD' * (M_h(MESH.internal_dof, MESH.internal_dof) * V_POD);
% t_assembly = toc(t_assembly);
% fprintf('done in %3.3f s', t_assembly);

%% Assemble stiffness matrix and K_star for the generalized alpha method
if DATA.DEIM_param_stiff
    theta_K = evaluate_ThetaFunctions_stiff(FOM_affine, param, DATA);
    K_in                   =  zeros(size(MESH.internal_dof,1), ...
                                    size(MESH.internal_dof,1));
    for i1 = 1 : FOM_affine.Qa
       K_in =  K_in + theta_K(1,i1).*FOM_affine.Aq{i1};
    end
    K_h =  zeros(MESH.numNodes*dim, MESH.numNodes*dim);
    K_h(MESH.internal_dof, MESH.internal_dof)  =  K_in;
    K   =  V_POD' * (K_in * V_POD);
    
        % per verificare la correttezza della decomposizione affine
%     K = SolidModel.compute_stiff_elastic(0, param); 
else
    K_h = SolidModel.compute_stiff_elastic();
    K   = V_POD' * (K_h(MESH.internal_dof, MESH.internal_dof) * V_POD);
end
K_star = Coef_Mass .* M + (1 - TimeAdvance.M_alpha_f) .* K;

%% Assemble external force vector F_ext
F_ext_FE = param(6) .* FOM_affine.F_ext;
F_ext    = V_POD' * F_ext_FE( MESH.internal_dof );

%% Initial Acceleration
F_ext_0  = F_ext;

F_in_0_FE  = K_h * u0;
F_in_0     = V_POD' * F_in_0_FE( MESH.internal_dof );

d2u0_N = M \ (F_ext_0 - F_in_0);

TimeAdvance.Initialize( V_POD'*u0(MESH.internal_dof), V_POD'*du0(MESH.internal_dof), d2u0_N);

LinSolver = LinearSolver( DATA.LinearSolver );

U_n_FE = u0(MESH.internal_dof);
U_n    = V_POD'*u0(MESH.internal_dof);

[~, ~, u_D]   =  CSM_ApplyEssentialBC([], [], MESH, DATA, t);

%% Time Loop
while ( t < tf )
    
    iter_time = tic;
    
    t       = t   + dt;
    k_t     = k_t + 1;
    
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n', t0, t, tf );

    t_force        = (1 - TimeAdvance.M_alpha_f)* t + ...
                     TimeAdvance.M_alpha_f * (t-dt);
    
    F_pre_FE       = FOM_affine.F_pressure .* param(3) .* sin(param(4).*t_force);
    F_Neu_FE       = FOM_affine.F_Neumann  .* param(6); 
    F_ext_FE       = FOM_affine.F_ext      .* param(5);
    
    F_pre          = V_POD' * F_pre_FE(MESH.internal_dof);
    F_Neu          = V_POD' * F_Neu_FE(MESH.internal_dof);
    F_ext          = V_POD' * F_ext_FE(MESH.internal_dof);
    
    Csi  = TimeAdvance.RhsContribute( );
    rhs  = - TimeAdvance.M_alpha_f * K * U_n + F_ext + M * Csi;
    
    rhs  = rhs + F_pre + F_Neu;
    
    % qualora non sia possibile una decompisizione affine delle condizioni
    % di Dirichlet
    % [~, ~, u_D]   =  CSM_ApplyEssentialBC([], [], MESH, DATA, t, 0);
    % U_k              = u(:,end);
    U_n_FE(MESH.Dirichlet_dof) = u_D .* param(7);
    
    Precon.Build( K_star );
    LinSolver.SetPreconditioner( Precon );
    
    %t_system = tic;
    U_n = LinSolver.Solve( K_star, rhs );
    %t_system = toc(t_system);
    %fprintf('done in %3.3f s\n', t_system); 
    
    U_n_FE(MESH.internal_dof)  = V_POD * U_n;
    
%     if sol_history
%         u = [u U_n];
%     else
%         u = U_n;
%     end
    
    
    %% Time advance and export monitored nodes displacements
    TimeAdvance.Update( U_n );
    
    if ~isempty(counter_mon_node_out)
        if counter_mon_node_out==0
            counter_node_done_out=counter_node_done_out+1;
            monitored_gdls.time_steps(counter_node_done_out) = t;
            monitored_gdls.displ(counter_node_done_out,:)= U_n_FE(dofID);
            counter_mon_node_out=n_mon_node_out;
        else
            counter_mon_node_out=counter_mon_node_out-1;
        end
    end
    
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
        [~, ~, ~] = mkdir(['/Users/matteotorzoni/Desktop/Corigliano/Mensola_matteo/Dati/Non_linear_fitting/Damage_Identification/mensola/damaged_',DATA.pb_ID_f,'/damaged_', DATA.pb_ID_f,'_', num2str(param(8))], 'graph_check');
        save(['/Users/matteotorzoni/Desktop/Corigliano/Mensola_matteo/Dati/Non_linear_fitting/Damage_Identification/mensola/damaged_', DATA.pb_ID_f, '/damaged_', DATA.pb_ID_f,'_', num2str(param(8)), '/graph_check/monitored_gdls_U_POD', num2str(index_param),'.mat']','-struct', 'monitored_gdls');
        save(['/Users/matteotorzoni/Desktop/Corigliano/Mensola_matteo/Dati/Non_linear_fitting/Damage_Identification/mensola/damaged_', DATA.pb_ID_f, '/damaged_', DATA.pb_ID_f,'_', num2str(param(8)), '/damage_level_1/monitored_gdls_U', num2str(i3),'.mat']','-struct', 'monitored_gdls');
    end
    %Questa non c'era
    [~, ~, ~] = mkdir(['/Users/matteotorzoni/Desktop/Corigliano/Mensola/Dati/Non_linear_fitting/Damage_Identification/mensola/damaged_', DATA.pb_ID_f, '/damaged_', DATA.pb_ID_f,'_', num2str(param(8)), '/damage_level_1']);
    %LA HO SPENTA IO
    %save(['/Users/matteotorzoni/Desktop/Corigliano/Mensola_matteo/Dati/Non_linear_fitting/Damage_Identification/mensola/damaged_', DATA.pb_ID_f, '/damaged_', DATA.pb_ID_f,'_', num2str(param(8)), '/damage_level_1/monitored_gdls_U', num2str(i3),'.mat']','-struct', 'monitored_gdls');
end

fprintf('\n************************************************************************* \n');

return
