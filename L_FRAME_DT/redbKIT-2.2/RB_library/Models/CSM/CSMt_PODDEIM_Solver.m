function [u, FE_SPACE, MESH, DATA, nodeID, dofID] = ...
    CSMt_PODDEIM_Solver(dim, elements, vertices,...
    boundaries, fem, data_file, param, vtk_filename, ...
    sol_history, ROM, node_positions, direction_nodeID, ...
    longer_analysis,n_mon_node_out, nodeID, dofID)
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

% 'special_save' allows to save sepatetly the time evolutions of the
    % basis (at maximum 4);
    % 'n_wanted_basis' specifies the number of basis that you may want to
    % save separately
special_save=0;
n_wanted_basis=2;

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
if longer_analysis
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
        dofID(i1)=nodeID(i1)+(direction_nodeID(i1))*MESH.numNodes;
    end
    monitored_gdls.time_steps=zeros(n_steps,1);
    monitored_gdls.displ=zeros(n_steps,length(nodeID));
    
    if special_save
        if n_wanted_basis >= 2
            monitored_gdls.displ1=zeros(n_steps,length(nodeID));
            monitored_gdls.displ2=zeros(n_steps,length(nodeID));
        end
        if n_wanted_basis >= 3
            monitored_gdls.displ3=zeros(n_steps,length(nodeID));
        end
        if n_wanted_basis >= 4
            monitored_gdls.displ4=zeros(n_steps,length(nodeID));
        end
    end
    
else
    dofID=[];
    nodeID=[];
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
F_in_0_FE  = SolidModel.compute_internal_forces( u0 );
F_in_0_FE  = F_in_0_FE( MESH.internal_dof );
F_in_0     = ROM.LeftProjection_int * F_in_0_FE( ROM.IDEIM_in );
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly)

d2u0_N = M \ (F_ext_0 - F_in_0);

TimeAdvance.Initialize( ROM.V'*u0(MESH.internal_dof), ROM.V'*du0(MESH.internal_dof), d2u0_N);%ROM.V'*d2u0(MESH.internal_dof) );

U_kN = ROM.V' * u(MESH.internal_dof,1);
U_n  = u0;
if special_save
    if n_wanted_basis >= 2
        u0_1 = zeros(size(u0));
        u0_2 = zeros(size(u0));
        u0_1(MESH.internal_dof) = V_POD(:,1)*U_kN(1);
        u0_2(MESH.internal_dof) = V_POD(:,2)*U_kN(2);
        u_1 = u0_1;
        u_2 = u0_2;
        U_n_1 = u0_1;
        U_n_2 = u0_2;
    end
    if n_wanted_basis >= 3
        u0_3 = zeros(size(u0));
        u0_3(MESH.internal_dof) = V_POD(:,3)*U_kN(3);
        u_3 = u0_3;
        U_n_3 = u0_3;
    end
    if n_wanted_basis >= 4
        u0_4 = zeros(size(u0));
        u0_4(MESH.internal_dof) = V_POD(:,4)*U_kN(4);
        u_4 = u0_4;
        U_n_4 = u0_4;
    end
end

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
    if special_save
        if n_wanted_basis >= 2
            dU_1             = zeros(MESH.numNodes*MESH.dim,1);
            dU_2             = zeros(MESH.numNodes*MESH.dim,1);
        end
        if n_wanted_basis >= 3
            dU_3             = zeros(MESH.numNodes*MESH.dim,1);
        end
        if n_wanted_basis >= 4
            dU_4             = zeros(MESH.numNodes*MESH.dim,1);
        end
    end
    U_k            = u(:,end);
    U_k(MESH.Dirichlet_dof) = u_D;
    if special_save
        if n_wanted_basis >= 2
            U_k_1            = u_1(:,end);
            U_k_2            = u_2(:,end);
            U_k_1(MESH.Dirichlet_dof) = 0;
            U_k_2(MESH.Dirichlet_dof) = 0;
        end
        if n_wanted_basis >= 3
            U_k_3            = u_3(:,end);
            U_k_3(MESH.Dirichlet_dof) = 0;
        end
        if n_wanted_basis >= 4
            U_k_4            = u_4(:,end);
            U_k_4(MESH.Dirichlet_dof) = 0;
        end
    end
    
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
    F_in_FE      = SolidModel.compute_internal_forces( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
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
        if special_save
            if n_wanted_basis >= 2
                dU_1(MESH.internal_dof) = V_POD(:,1)*dU_N(1);
                dU_2(MESH.internal_dof) = V_POD(:,2)*dU_N(2);
            end
            if n_wanted_basis >= 3
                dU_3(MESH.internal_dof) = V_POD(:,3)*dU_N(3);
            end
            if n_wanted_basis >= 4
                dU_4(MESH.internal_dof) = V_POD(:,4)*dU_N(4);
            end
        end
        fprintf('\n        time to solve the linear system in %3.3f s \n', LinSolver.GetSolveTime());
        
        % update solution
        U_k        = U_k + dU;
        if special_save
            if n_wanted_basis >= 2
                U_k_1 = U_k_1 + dU_1;
                U_k_2 = U_k_2 + dU_2;
            end
            if n_wanted_basis >= 3
                U_k_3 = U_k_3 + dU_3;
            end
            if n_wanted_basis >= 4
                U_k_4 = U_k_4 + dU_4;
            end
        end
        U_kN       = U_kN + dU_N;
        incrNorm   = norm(dU)/norm(U_k);
        
        % Assemble matrix and right-hand side
        fprintf('\n   -- Assembling internal forces... ');
        t_assembly = tic;
        F_in_FE      = SolidModel.compute_internal_forces( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
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
        if special_save
            if n_wanted_basis >= 2
                u_1 = U_k_1;
                u_2 = U_k_2;
            end
            if n_wanted_basis >= 3
                u_3 = U_k_3;
            end
            if n_wanted_basis >= 4
                u_4 = U_k_4;
            end
        end
    end
    
    if ~isempty(counter_mon_node_out)
        if counter_mon_node_out==0
            counter_node_done_out=counter_node_done_out+1;
            monitored_gdls.time_steps(counter_node_done_out) = t;
            monitored_gdls.displ(counter_node_done_out,:)= U_n(dofID);
            if special_save
                if n_wanted_basis >= 2
                    monitored_gdls.displ1(counter_node_done_out,:)= U_n_1(dofID);
                    monitored_gdls.displ2(counter_node_done_out,:)= U_n_2(dofID);
                    check = U_n_1(dofID) + U_n_2(dofID);
                end
                if n_wanted_basis >= 3
                    monitored_gdls.displ3(counter_node_done_out,:)= U_n_3(dofID);
                    check = check + U_n_3(dofID);
                end
                if n_wanted_basis >= 4
                    monitored_gdls.displ4(counter_node_done_out,:)= U_n_4(dofID);
                    check = check + U_n_4(dofID);
                end
                check = check - U_n(dofID);
                if isempty(check)
                    disp('*** warning - the decomposed solution is wrong')
                    return
                end
            end
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
    monitored_gdls.time_steps=monitored_gdls.time_steps(1:number_entries);
    monitored_gdls.displ=monitored_gdls.displ(1:number_entries,:);
    if special_save
        if n_wanted_basis >= 2
            monitored_gdls.displ1=monitored_gdls.displ1(1:number_entries,:);
            monitored_gdls.displ2=monitored_gdls.displ2(1:number_entries,:);
        end
        if n_wanted_basis >= 3
            monitored_gdls.displ3=monitored_gdls.displ3(1:number_entries,:); %#ok<*STRNU>
        end
        if n_wanted_basis >= 4
            monitored_gdls.displ4=monitored_gdls.displ4(1:number_entries,:); %#ok<*STRNU>
        end
    end
    save('C:\Users\Rosafalco\Documents\Polimi\3_Dottorato\Codici\Matlab\Luca\graphkit\DFT\monitored_gdls_U_POD_POD.mat','-struct', 'monitored_gdls');
    save('C:\Users\Rosafalco\Documents\Polimi\3_Dottorato\Codici\Matlab\Luca\graphkit\General\monitored_gdls_U_POD_POD.mat','-struct', 'monitored_gdls');
end

fprintf('\n************************************************************************* \n');

return
