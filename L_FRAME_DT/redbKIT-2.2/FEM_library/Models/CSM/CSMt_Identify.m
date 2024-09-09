function [dofID, monitored_gdls] = CSMt_Istanze_Damping(MESH, FE_SPACE, DATA, param, node_positions, direction_nodeID, ...
        dofID, Labeling, figura, str_out)

%CSMT_SOLVER Elasto-Dynamic Structural

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

TimeAdvance = GeneralizedAlpha_TimeAdvance_Damped( DATA.time.beta, DATA.time.gamma, DATA.time.alpha_m, DATA.time.alpha_f, dt );

Coef_Mass = TimeAdvance.MassCoefficient( );
Coef_Damp = TimeAdvance.DampingCoefficient( );
Coef_Stif = TimeAdvance.StiffnessCoefficient( );

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

PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Preconditioner.type, DATA);

%% Creo modello
SolidModel = CSM_Assembler_matteo( MESH, DATA, FE_SPACE );

%% Assemble stiffness matrix and K_star for the generalized alpha method
K = SolidModel.compute_stiff_elastic(0, param);
M   = SolidModel.compute_mass_spatial();
C = sparse(zeros(size(K)));
C(MESH.internal_dof,MESH.internal_dof) = SolidModel.compute_damping( K(MESH.internal_dof,MESH.internal_dof), M(MESH.internal_dof,MESH.internal_dof), C(MESH.internal_dof,MESH.internal_dof));
K_star = Coef_Mass .* M + Coef_Damp .* C + Coef_Stif .* K;

LinSolver = LinearSolver( DATA.LinearSolver );
Precon.Build( K_star(MESH.internal_dof,MESH.internal_dof) );
LinSolver.SetPreconditioner( Precon );

%% Assemble Neumann and pressure force vector
Forzanti       = Valuta_forzanti(FE_SPACE, MESH, DATA, param, 0);

%% Assemble external force vector F_ext
%F_inerzia    = SolidModel.compute_inertia_forces(0, param);

%% Initial Acceleration
F_ext_0  = Forzanti.F_Neumann .* 0;% + F_inerzia .* ACC(1+round(t/dt));
F_in_0   =  K * u0;

d2u0=zeros(MESH.numNodes*dim,1);
d2u0(MESH.internal_dof) = M(MESH.internal_dof,MESH.internal_dof) \ (F_ext_0(MESH.internal_dof) - F_in_0(MESH.internal_dof));

TimeAdvance.Initialize( u0, du0, d2u0 );

U_n = u0;

[~, ~, u_D]   =  CSM_ApplyEssentialBC([], [], MESH, DATA, t);

%% Time Loop
while ( (t+0.00001) < tf )

    t       = t   + dt;
    
    %fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n', t0, t, tf );

    %t_force        = (1 - TimeAdvance.M_alpha_f)* t + TimeAdvance.M_alpha_f * (t-dt);
   
    F_pre          = Forzanti.F_pressure .* param(3) .* sin(param(4).*2*3.14159.*t);%t_force
    
    Csi  = TimeAdvance.RhsContribute_M( );
    Phi  = TimeAdvance.RhsContribute_C( );
    Csa  = TimeAdvance.RhsContribute_K( );
    
    rhs  = F_pre + (M * Csi) + (C * Phi) + (K * Csa);
    
    U_n(MESH.Dirichlet_dof) = u_D;
    
    U_n(MESH.internal_dof) = ...
    LinSolver.Solve( K_star(MESH.internal_dof,MESH.internal_dof), ...
                     rhs(MESH.internal_dof));

    %% Time advance and export monitored nodes displacements
    TimeAdvance.Update( U_n );
    
    %U2_n = TimeAdvance.M_d2U;
    if t >= t_start
        %monitored_gdls.acc(round(1+(t-t_start)/dt),:)= U2_n(dofID);
        monitored_gdls.displ(round(1+(t-t_start)/dt),:)= U_n(dofID);
    end
    if figura
        WW(jj)=U_n(3242);
        %U2_n=TimeAdvance.M_d2U;
        %WW2(jj)=U2_n(3161);
        %WW(jj)=U_n(3161);
        jj=jj+1;
    end
    
end

%% Plot
if figura
    x=1:1:(round(tf/dt));
    figure
    plot(x.*dt,WW)
    xlim([0 tf])
%     figure
%     plot(x.*dt,WW2)
%     xlim([0 tf])
end
    
%% Save and Labeling

monitored_gdls.Amplitude = Labeling(1);
monitored_gdls.Frequency = Labeling(2);
monitored_gdls.Damage_class = Labeling(3);
monitored_gdls.Damage_level = Labeling(4);

for i=1:size(monitored_gdls.displ,2)
    %dlmwrite([str_out,'/U2_concat_',num2str(i),'.csv'],monitored_gdls.acc(:,i),'-append')
    dlmwrite([str_out,'/U_concat_',num2str(i),'.csv'],monitored_gdls.displ(:,i),'-append')
end
dlmwrite([str_out,'/Amplitude.csv'],monitored_gdls.Amplitude,'-append')
dlmwrite([str_out,'/Frequency.csv'],monitored_gdls.Frequency,'-append')
dlmwrite([str_out,'/Damage_class.csv'],monitored_gdls.Damage_class,'-append')
dlmwrite([str_out,'/Damage_level.csv'],monitored_gdls.Damage_level,'-append')

return
