function [dofID] = CSMt_POD_Solver_Spatial(MESH, FE_SPACE, ...
        DATA, param, node_positions, direction_nodeID, ...
        dofID, V_POD, pb_numb, white_noise_corruption, SNR, Labeling, M, figura, path)

%CSMT_POD_SOLVER Elasto-Dynamic Structural

dim=FE_SPACE.dim;
DATA.param = param;

if figura
    jj=1;
end

%% Gather Time Setting
t0        = DATA.time.t0;
dt        = DATA.time.dt;
tf        = DATA.time.tf_long;
t         = DATA.time.t0;
t_start   = DATA.time.t_start;

TimeAdvance = GeneralizedAlpha_TimeAdvance( DATA.time.beta, DATA.time.gamma, DATA.time.alpha_m, DATA.time.alpha_f, dt );
Coef_Mass = TimeAdvance.MassCoefficient( );

%% Determine the dofs to be saved
if  isempty(dofID)    
    % estimate the number of required steps
    n_steps=round(1+(tf-t_start)/dt);   
    % find the ID of the monitored nodes
    nodeID=ones(1,size(node_positions,2));
    for j=1:size(node_positions,2)
        translated_vertices = MESH.vertices - node_positions(:,j);
        norm_translation=100.*ones(1,size(translated_vertices,2));
        for i2 = 1:length(norm_translation)
            if MESH.dim == 2
                norm_translation(i2) = (translated_vertices(1,i2))^2+(translated_vertices(2,i2))^2;
            elseif MESH.dim == 3
                norm_translation(i2) = (translated_vertices(1,i2))^2+(translated_vertices(2,i2))^2+(translated_vertices(3,i2))^2;
            end
        end
        [~,nodeID(j)]=min(norm_translation);
    end

    % find the ID of the monitored dofs
    dofID=zeros(1,size(node_positions,2));
    for j=1:length(nodeID)
        dofID(j) = nodeID(j)+(direction_nodeID(j))*MESH.numNodes;
    end
    
    monitored_gdls.displ=zeros(n_steps,length(dofID));
    monitored_gdls.acc=zeros(n_steps,length(dofID));
else
    n_steps=round(1+(tf-t_start)/dt); 
    monitored_gdls.displ=zeros(n_steps,length(dofID));
    monitored_gdls.acc=zeros(n_steps,length(dofID));
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

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf(' * Number of Dofs      = %d \n',size(V_POD, 2));
fprintf(' * Number of timesteps =  %d\n', (tf-t0)/dt);
fprintf('-------------------------------------------\n');

PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Preconditioner.type, DATA);

%% Creo modello
SolidModel = CSM_Assembler_matteo( MESH, DATA, FE_SPACE );

%% Assemble stiffness matrix and K_star for the generalized alpha method
K_h      = SolidModel.compute_stiff_elastic(0, param);
K   =  V_POD' * (K_h(MESH.internal_dof,MESH.internal_dof) * V_POD);
K_star = Coef_Mass .* M + (1 - TimeAdvance.M_alpha_f) .* K;

LinSolver = LinearSolver( DATA.LinearSolver );
Precon.Build( K_star );
LinSolver.SetPreconditioner( Precon );

%% Assemble Neumann and pressure force vector
Forzanti       = Valuta_forzanti(FE_SPACE, MESH, DATA, param, 0);

%% Initial Acceleration
F_ext_0_FE  = Forzanti.F_Neumann .* 0;
F_ext_0    = V_POD' * F_ext_0_FE( MESH.internal_dof );

F_in_0_FE  = K_h * u0;
F_in_0     = V_POD' * F_in_0_FE( MESH.internal_dof );

d2u0 = M \ (F_ext_0 - F_in_0);

TimeAdvance.Initialize( V_POD'*u0(MESH.internal_dof), V_POD'*du0(MESH.internal_dof), d2u0);

U_n_FE = u0;
U_n    = V_POD'*u0(MESH.internal_dof);

[~, ~, u_D]   =  CSM_ApplyEssentialBC([], [], MESH, DATA, t);

%% Time Loop
while ( (t+0.00001) < tf )
    
    t       = t   + dt;
    
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n', t0, t, tf );

    t_force = (1 - TimeAdvance.M_alpha_f)* t + TimeAdvance.M_alpha_f * (t-dt);
                 
    F_pre_FE       = Forzanti.F_pressure .* param(3) .* sin(param(4).*2*3.14159.*t_force); 
    F_pre          = V_POD' * F_pre_FE(MESH.internal_dof);
    
    Csi  = TimeAdvance.RhsContribute( );
    rhs  = - TimeAdvance.M_alpha_f * K * U_n + M * Csi + F_pre;
    
    U_n_FE(MESH.Dirichlet_dof) = u_D;
    
    U_n = LinSolver.Solve( K_star, rhs );
    
    U_n_FE(MESH.internal_dof)  = V_POD * U_n;
    
    %% Time advance and export monitored nodes displacements
    TimeAdvance.Update( U_n );
    
    acc_ridotte = TimeAdvance.M_d2U;
    U_2_FE(MESH.internal_dof) = V_POD * acc_ridotte;
    if t >= t_start
        monitored_gdls.acc(round(1+(t-t_start)/dt),:)= U_2_FE(dofID);
        monitored_gdls.displ(round(1+(t-t_start)/dt),:)= U_n_FE(dofID);
    end
    if figura
        WW2(jj)=U_2_FE(3161);
        WW(jj)=U_n_FE(3161);
        jj=jj+1;
    end
end

%% Plot
if figura
    x=1:1:(round(tf/dt));
    figure
    hold on
    plot(x.*dt,WW)
    xlim([0 tf])
    box on
    figure
    hold on
    plot(x.*dt,WW2)
    xlim([0 tf])
    box on
end

%% Labeling and Save
if Labeling(1,1)<3.7
    monitored_gdls.label_x = 0.3/2+Labeling(1,1);
    monitored_gdls.label_y = 0.3/2;
elseif Labeling(1,1)>=3.7
    monitored_gdls.label_x = 4-0.3/2;
    monitored_gdls.label_y = 0.3/2+Labeling(1,1)-3.7;
end

monitored_gdls.label_d=Labeling(2,1);

%Se non serve senza rumore commentare
% str_out   = [path,'/Dati/istantrain_',pb_numb];
% for i=1:size(monitored_gdls.acc,2)
%     dlmwrite([str_out,'/U_concat_',num2str(i),'.csv'],monitored_gdls.displ(:,i),'-append')
%     dlmwrite([str_out,'/U2_concat_',num2str(i),'.csv'],monitored_gdls.acc(:,i),'-append')
% end
% dlmwrite([str_out,'/label_x.csv'],monitored_gdls.label_x,'-append')
% dlmwrite([str_out,'/label_y.csv'],monitored_gdls.label_y,'-append')
% dlmwrite([str_out,'/label_d.csv'],monitored_gdls.label_d,'-append')
   

%% wgn corruption
if white_noise_corruption
    for sensore=1:size(dofID,2)
        % calcoli lo scarto quadratico medio
        rms_signal = rms(abs(monitored_gdls.displ(:,sensore)));
        % calcoli la deviazione standard necessaria per avere rumore bianco
        dev_std    = rms_signal / sqrt(SNR);   
        % funzione che aggiunge rumore bianco
        [monitored_gdls.displ(:,sensore)] = white_gauss_noise(monitored_gdls.displ(:,sensore), dev_std);
        % calcoli lo scarto quadratico medio
        rms_signal = rms(abs(monitored_gdls.acc(:,sensore)));  
        % calcoli la deviazione standard necessaria per avere rumore bianco
        dev_std    = rms_signal / sqrt(SNR);      
        % funzione che aggiunge rumore bianco
        [monitored_gdls.acc(:,sensore)] = white_gauss_noise(monitored_gdls.acc(:,sensore), dev_std);
    end    

    str_out   = [path,'/Dati/istantrain_',pb_numb,'_SNR_',num2str(SNR)];
    for i=1:size(monitored_gdls.acc,2)
        dlmwrite([str_out,'/U_concat_',num2str(i),'.csv'],monitored_gdls.displ(:,i),'-append')
        dlmwrite([str_out,'/U2_concat_',num2str(i),'.csv'],monitored_gdls.acc(:,i),'-append')
    end
    dlmwrite([str_out,'/label_x.csv'],monitored_gdls.label_x,'-append')
    dlmwrite([str_out,'/label_y.csv'],monitored_gdls.label_y,'-append')
    dlmwrite([str_out,'/label_d.csv'],monitored_gdls.label_d,'-append')
end
fprintf('\n************************************************************************* \n');

return
