function [FE_SPACE, MESH, DATA] = ...
                         CSMt_Solver_El_BASI(MESH, FE_SPACE,...
                         DATA, param, Training_Options, figura)

%CSMT_SOLVER Elasto-Dynamic Structural Solver

dim=FE_SPACE.dim;
DATA.param=param;

if figura
    jj=1;
end

%% Gather Time Setting
t0        = DATA.time.t0;
dt        = DATA.time.dt;
tf        = DATA.time.tf_short;
t         = DATA.time.t0;

TimeAdvance = GeneralizedAlpha_TimeAdvance_Damped( DATA.time.beta, DATA.time.gamma, DATA.time.alpha_m, DATA.time.alpha_f, dt );

Coef_Mass = TimeAdvance.MassCoefficient( );
Coef_Damp = TimeAdvance.DampingCoefficient( );
Coef_Stif = TimeAdvance.StiffnessCoefficient( );

%% Variable initialization

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
fprintf(' * Number of Dofs      = %d \n',length(MESH.internal_dof));
fprintf(' * Number of timesteps =  %d\n', (tf-t0)/dt);
fprintf('-------------------------------------------\n');

DispSnap_h5  = HDF5_DenseMultiCVector(Training_Options.h5_filename, Training_Options.h5_section, length(MESH.internal_dof));
DispSnap_h5.append( u0(MESH.internal_dof) );

PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Preconditioner.type, DATA);

%% carico dati per la temperatura già risolta in casi random
SolidModel = CSM_Assembler_matteo( MESH, DATA, FE_SPACE );

%% Assemble stiffness matrix and K_star for the generalized alpha method 
K = SolidModel.compute_stiff_elastic(0, param);
M = SolidModel.compute_mass_spatial();
C = sparse(zeros(size(K)));
C(MESH.internal_dof,MESH.internal_dof) = SolidModel.compute_damping( K(MESH.internal_dof,MESH.internal_dof), M(MESH.internal_dof,MESH.internal_dof), C(MESH.internal_dof,MESH.internal_dof));
K_star = Coef_Mass .* M + Coef_Damp .* C + Coef_Stif .* K;

Precon.Build( K_star(MESH.internal_dof,MESH.internal_dof) );
LinSolver = LinearSolver( DATA.LinearSolver );
LinSolver.SetPreconditioner( Precon );

%% Assemble Neumann and pressure force vector
Forzanti       = Valuta_forzanti(FE_SPACE, MESH, DATA, param, 0);

%% Assemble external force vector F_ext
%F_inerzia    = SolidModel.compute_inertia_forces(0, param);

%% Initial Acceleration
F_ext_0  = Forzanti.F_Neumann .* 0;
F_in_0   =  K * u0;

%considera che l'imposizione di forze iniziali sia equivalente
%all'imposizione di un'accelerazione inziale
d2u0=zeros(MESH.numNodes*dim,1);
d2u0(MESH.internal_dof) = M(MESH.internal_dof,MESH.internal_dof) \ (F_ext_0(MESH.internal_dof) - F_in_0(MESH.internal_dof));

TimeAdvance.Initialize( u0, du0, d2u0 );

U_n = u0;

[~, ~, u_D]   =  CSM_ApplyEssentialBC([], [], MESH, DATA, t);

%% Time Loop
while ( (t+0.00001) < tf )
    
    t       = t   + dt;
    
    %fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n', t0, t, tf );
    
    %t_force = (1 - TimeAdvance.M_alpha_f)* t + TimeAdvance.M_alpha_f * (t-dt);
    
    F_pre          = Forzanti.F_pressure .* param(3) .* sin(param(4).*2*3.14159.*t);%t_force
    
    Csi  = TimeAdvance.RhsContribute_M( );
    Phi  = TimeAdvance.RhsContribute_C( );
    Csa  = TimeAdvance.RhsContribute_K( );
    
    rhs  = F_pre + (M * Csi) + (C * Phi) + (K * Csa);
    
    U_n(MESH.Dirichlet_dof) = u_D;
    
    U_n(MESH.internal_dof) = ...
    LinSolver.Solve( K_star(MESH.internal_dof,MESH.internal_dof), ...
                     rhs(MESH.internal_dof));
                 
    DispSnap_h5.append( U_n(MESH.internal_dof) );
    
    %% Time advance & monitored gdls   
    TimeAdvance.Update( U_n );
    
    if figura
        %U2_n=TimeAdvance.M_d2U;
        WW(jj)=U_n(3242);
        %WW(jj)=U_n(3161);
        %WW2(jj)=U2_n(3161);
        jj=jj+1;
    end
    
end

%% Plot monitored gdls

    if figura
        x=1:1:(round(tf/dt));
        figure
        hold on
        plot(x.*dt,WW)
        xlim([0 tf])
        box on
%         figure
%         hold on
%         plot(x.*dt,WW2)
%         xlim([0 tf])
%         box on
    end

fprintf('\n************************************************************************* \n');
return
