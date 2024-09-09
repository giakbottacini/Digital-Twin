clc
clear all
%%
ESEMPIO = 'L_FRAME_DT';
cd D:\Users\Matteo\Corigliano\L_FRAME_DT
path = ['D:\Users\Matteo\Corigliano\',ESEMPIO];

%% DATA
fem      =  'P1'; %CST
dim      = 3; %2Dproblem
mesh_file       = 'L_FRAME';
[vertices, boundaries, elements] = msh_to_Mmesh(['mesh\',mesh_file], dim);

% per caricare il modello già salvato:
%   -(load_FOM = true)
% per calcolare le autofrequenze:
%   -(calcolo_autofrequenze = true)
% per plottare i modi propri in vtk:
%   -(stampa_modi = true)
% per caricare il sampling dello spazio parametrico salvato:
%   -(load_lhs = true)
% per allenare modello ridotto per problema meccanico:
%   -(load_ROM_displ = 0)
% per generare istanze per allenare la rete:
%   -(train_instance = true)
% per corrompere con del rumore bianco le istanze di train:
%   -(white_noise_corruption = true)
% per generare istanze per testare la rete (FOM):
%   -(test_instance = true)
% per plottare figura di prova:
%   -(figura = true)

load_FOM               = 1;
calcolo_autofrequenze  = true;
    stampa_autovettori = false;
    num_autovettori    = 10;
load_lhs               = 1;
load_ROM_displ         = 1;
train_instance         = false;
test_instance          = false;
figura                 = false;

identify_genetic       = false;

pb_ID_f         = '7_1'; % Distinguo tra i vari modelli ridotto
pb_numb         = '7_1'; % Distinguo per la creazione istanze
pb_identifiy    = '7_1';

data_file       = 'CSMGL_data_multi';

FOM_n     = 7; % File ID for the training parameters

% coordinate dei nodi di cui si richiede l'uscita: gdl monitorati 
node_positions         = [1, 2, 3, 3.85, 3.85, 3.85, 3.85, 3.85;...
                          0.15, 0.15, 0.15, 0.15, 1, 2, 3, 4;...
                          0, 0, 0, 0, 0, 0, 0, 0];
node_directions        = [2, 2, 2, 2, 2, 2, 2, 2];

dofID  = [];
% Set quad_order
if dim == 2
    quad_order       = 4; %Ordine di integrazione (dovrebbe funzionare anche con 2?)
elseif dim == 3
    quad_order       = 5;
end

% Read problem data
DATA   = CSM_read_DataFile(data_file, dim); %prende parametri per il FOM

%% FE_space and Mesh
if ~load_FOM
[ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA, 'CSM' );
[ FE_SPACE ] = buildFESpace( MESH, fem, dim, quad_order );
        [~, ~, ~] = mkdir([path,'\Dati\'],['damaged_',pb_ID_f]);
        save([path,'\Dati\damaged_', pb_ID_f,'\MESH.mat']','-struct', 'MESH');
        save([path,'\Dati\damaged_', pb_ID_f,'\FE_SPACE.mat']','-struct', 'FE_SPACE');
else
    str_FOM_in = [path,'\Dati\damaged_', pb_ID_f,];
    MESH              = load([str_FOM_in,'\MESH.mat']);
    FE_SPACE          = load([str_FOM_in,'\FE_SPACE.mat']);
end

%% Calcolo autofrequenze
if calcolo_autofrequenze

    Default_Parameters  = [30e+09;  ...
                            0.2;    ...
                            0;      ...
                            0;
                            0;
                            0];
    
    Default_Parameters = CSM_correct_param(DATA,dim,Default_Parameters); 
                            
    SolidModel = CSM_Assembler_matteo( MESH, DATA, FE_SPACE );
    K      = SolidModel.compute_stiff_elastic(0, Default_Parameters);
    K_in   =  K(MESH.internal_dof,MESH.internal_dof);
    M_h   = SolidModel.compute_mass_spatial();
    M_in = M_h(MESH.internal_dof,MESH.internal_dof);
    [eigenmod_flex, eigval] = eigs(K_in,M_in,num_autovettori,'smallestabs');
    eigval = diag(eigval);
    omega  = sqrt(eigval);
    eigen = omega / (2*3.14159);
    for i=1:size(eigen,1)
        autoval(i)=real(eigen(i));
    end

    %Plot dei autovettori per il problema meccanico
    if stampa_autovettori
        forma=zeros(size(vertices,2)*dim,num_autovettori);
        for i=1:num_autovettori
            forma(MESH.internal_dof,i)=eigenmod_flex(:,i);
            CSM_export_solution(MESH.dim, forma(:,i), MESH.vertices, MESH.elements, MESH.numNodes, ['Forma modale ',num2str(i)]);
        end
    return
    end
end

%% Parameters
if ~load_lhs
    lhs=lhsdesign(DATA.how_many_snaps,4)'; 
    Training_Parameters_FOM=zeros(6,DATA.how_many_snaps);
    Training_Parameters_FOM(1,:)=ones(1,DATA.how_many_snaps).*30e09;
    Training_Parameters_FOM(2,:)=ones(1,DATA.how_many_snaps).*0.2;
    Training_Parameters_FOM(3,:)=ones(1,DATA.how_many_snaps).*40000+lhs(1,:).*(80000-40000);
    Training_Parameters_FOM(4,:)=ones(1,DATA.how_many_snaps).*10+lhs(2,:).*(60-10);%.*5+lhs(2,:).*(90-5)
    Training_Parameters_FOM(5,:)=floor(lhs(3,:).*(8));
    Training_Parameters_FOM(6,:)=ones(1,DATA.how_many_snaps).*0.3+lhs(4,:).*(0.8-0.23);%.*0.1+floor(lhs(4,:).*(8)).*0.1;
    [~, ~, ~] = mkdir([path,'\redbKIT-2.2\'],'train_param');
    csvwrite([path,'\redbKIT-2.2\train_param\param_FOM_',num2str(FOM_n),'.csv'],Training_Parameters_FOM)
end

%% ROM generation/loading displacements
if ~load_ROM_displ % GENERO BASI RIDOTTE PER GLI SPOSTAMENTI
    
    [~,~,~] = mkdir(['Snapshots\',num2str(pb_ID_f)]);
    
    % Parameter loading 
    Training_Parameters_FOM = csvread([path,'\redbKIT-2.2\train_param\param_FOM_',num2str(FOM_n),'.csv']);
    % Parameter elaboration
    Training_Parameters = CSM_correct_param(DATA,dim,Training_Parameters_FOM); %passo in costanti di lamè + correzione per planestress
    
    for i = 1 : DATA.how_many_snaps
        delete(['Snapshots\',num2str(pb_ID_f),'\DisplacementSnapshotsMP',num2str(i),'.h5'])
    end
    
    % inizializzazione delle basi dati
    h5_filename_Sol       = ['Snapshots\',num2str(pb_ID_f),'\DisplacementSnapshotsMP'];
    OfflineTraining.Solution.h5_section        = 'Displacement';
    OfflineTraining.Solution.SamplingFrequency = 1;
    V = [];
    
    for i = 1 : DATA.how_many_snaps
        OfflineTraining.Solution.h5_filename = [h5_filename_Sol, num2str(i), '.h5'];
        [~ , ~, ~] = CSMt_Solver_El_BASI(MESH, FE_SPACE,...
         DATA, Training_Parameters(:,i), OfflineTraining.Solution, figura);
    end
    for i = 1 : DATA.how_many_snaps
        DispSnap_HDF5  = HDF5_DenseMultiCVector([h5_filename_Sol, num2str(i), '.h5'], OfflineTraining.Solution.h5_section);
        S_u            = DispSnap_HDF5.readValues();
        S_size         = size(S_u        ,2);
        [V_local, ~, Sigma_loc] = VPOD_basis_computation(S_u,         [], DATA.tol_POD_U_local, 1);
        S_size         = size([V V_local],2);
        [V      , ~, Sigma    ] = VPOD_basis_computation([V V_local], [], DATA.tol_POD_U,       1);
    end    
    ROM_displ.V = V;
    save([path,'\Dati\damaged_', pb_ID_f, '\ROM_displ.mat']','-struct', 'ROM_displ');
end  
%% Genero istanze train
if train_instance == true
    
    str_ROM_in = [path,'\Dati\damaged_', pb_ID_f];
    ROM_displ = load([str_ROM_in,'\ROM_displ.mat']);  
    
    lhs=lhsdesign(DATA.how_many_train,4)'; 
    Param_case=zeros(4,DATA.how_many_train);
    Param_case(1,:)=ones(1,DATA.how_many_train).*40000+lhs(1,:).*(80000-40000);
    Param_case(2,:)=ones(1,DATA.how_many_train).*10+lhs(2,:).*(60-10);%.*5+lhs(2,:).*(90-5)
    Param_case(3,:)=floor(lhs(3,:).*(8));
    Param_case(4,:)=ones(1,DATA.how_many_train).*0.3+lhs(4,:).*(0.8-0.3);%.*0.1+floor(lhs(4,:).*(8)).*0.1;
    
    str_out   = [path,'\Dati\istantrain_',pb_numb];
    [~, ~, ~] = mkdir(str_out);
    
    for i1 = 1 : DATA.how_many_train
        
        Instance_Parameters  = [30e+09;                             ... 
                                0.2;                                ...
                                Param_case(1,i1);                   ...
                                Param_case(2,i1);                   ...
                                Param_case(3,i1);                   ...
                                Param_case(4,i1)];
            
        Istance_Parameters = CSM_correct_param(DATA,dim,Instance_Parameters); 
        
        [dofID] = CSMt_Istanze_Damping_ROM(MESH, FE_SPACE, DATA, Istance_Parameters, node_positions, node_directions, ...
        dofID, ROM_displ.V, Param_case(:,i1), figura, str_out);

    end
end

%% Genero istanze test
if test_instance == true
    lhs=lhsdesign(DATA.how_many_test,4)'; 
    Param_case=zeros(4,DATA.how_many_test);
    Param_case(1,:)=ones(1,DATA.how_many_test).*40000+lhs(1,:).*(80000-40000);
    Param_case(2,:)=ones(1,DATA.how_many_test).*10+lhs(2,:).*(60-10);%.*5+lhs(2,:).*(90-5)
    Param_case(3,:)=floor(lhs(3,:).*(8));
    Param_case(4,:)=ones(1,DATA.how_many_test).*0.3+lhs(4,:).*(0.8-0.3);%floor(lhs(4,:).*(8)).*0.1;
    
    str_out   = [path,'\Dati\istantest_',pb_numb];
    [~, ~, ~] = mkdir(str_out);
   
    for i1 = 1 : DATA.how_many_test

        Instance_Parameters  = [30e+09;                             ... 
                                0.2;                                ...
                                Param_case(1,i1);                   ...
                                Param_case(2,i1);                   ...
                                Param_case(3,i1);                   ...
                                Param_case(4,i1)];

        Istance_Parameters = CSM_correct_param(DATA,dim,Instance_Parameters); 

        [dofID] = CSMt_Istanze_Damping(MESH, FE_SPACE, DATA, Istance_Parameters, node_positions, node_directions, ...
        dofID, Param_case(:,i1), figura, str_out);
    end
end

%% Identify
if identify_genetic == true
    
    str_ROM_in = [path,'\Dati\damaged_', pb_ID_f];
    ROM_displ = load([str_ROM_in,'\ROM_displ.mat']);  
    str_out   = [path,'\Dati\istanidentify_',pb_identifiy];
    [~, ~, ~] = mkdir(str_out);
    
    lhs=lhsdesign(DATA.how_many_identify,4)'; 
    Param_case=zeros(4,DATA.how_many_identify);
    Param_case(1,:)=ones(1,DATA.how_many_identify).*40000+lhs(1,:).*(80000-40000);
    Param_case(2,:)=ones(1,DATA.how_many_identify).*10+lhs(2,:).*(60-10);
    Param_case(3,:)=floor(lhs(3,:).*(8));
    Param_case(4,:)=ones(1,DATA.how_many_identify).*0.3+lhs(4,:).*(0.8-0.3);
    
    min_x = 10;
    max_x = 60;
    min_y = 0.3;
    max_y = 0.8;
    
    % To treat the discrete parameter controling the damage zone we perform
    % a brute force optimization. To deal with the load amplitude parameter
    % we simply consder the structure as linear and hence we assume to
    % already know the associated value of the parameter; this wuold be
    % equivalent to normalize the observed and the output signals. We
    % consider signals not corrupted with noise, to facilitate the comparison.
    % The main bottolneck is given by the computing time reqired to update
    % the reduced order model; on the other hand NNs are very effective
    % during the online phase; the same is also true, for instance, with
    % decision trees, but these latter require to simulate offline all the
    % possible configurations, which is not feasible in general; moreover 
    % this latter would probably require the use of some damage sensitive 
    % feature in place of the raw vibration recordings. On the other hand, 
    % the problem with the use of damage-sensitive features is that a model 
    % updating strategy directly exploiting the ROM would fail to update the 
    % parameters ruling the loading condition if these do not affect the 
    % adopted features; a possible solutin would be to update a surrogate 
    % model instead, which map the parameters of interest onto sensitive-features.
    
    % The optimization takes around 1 minute and 15 seconds / 100iterations. 
    % we run optimizations with 200 iterations, so that considering the 8
    % possible subdomains we end up with 1200 iterations required per
    % instance; this results into a computing time of 20 minutes for each
    % instance. 1600/100*75/60
    
    for i1 = 1 : DATA.how_many_identify

        Instance_Parameters  = [30e+09;                             ... 
                                0.2;                                ...
                                Param_case(1,i1);                   ...
                                Param_case(2,i1);                   ...
                                Param_case(3,i1);                   ...
                                Param_case(4,i1)];

        Istance_Parameters = CSM_correct_param(DATA,dim,Instance_Parameters); 
        
        [dofID, observed] = CSMt_Identify(MESH, FE_SPACE, DATA, Istance_Parameters, node_positions, node_directions, ...
        dofID, Param_case(:,i1), figura, str_out);
        
        for i2 = 1 : 8 %cycle on the number of subdomains
            iteration = 1;
            for i3 = 1:4 %cycle on the initial population
                
                %sample parameters
                x(i3) = (max_x - min_x)*rand(1,1) + min_x;
                y(i3) = (max_y - min_y)*rand(1,1) + min_y;
                
                Output_Parameters(i3,:)  = Istance_Parameters;
                Output_Parameters(i3,4)  = x(i3);
                Output_Parameters(i3,5)  = i2-1;
                Output_Parameters(i3,6)  = y(i3);
                
                %eval ROM
                [dofID,output] = CSMt_Identify_ROM(MESH, FE_SPACE, DATA, Output_Parameters(i3,:), node_positions, node_directions, ...
                dofID, ROM_displ.V, Param_case(:,i1), figura, str_out);
    
                %fitness evaluation
                fit(i3) = 0.5 * norm(observed.displ - output.displ,'fro')^2 ; %+ np.linalg.norm(x)
                    
            end
            
            while iteration < 200
                
                fprintf('\n************************************************************************* \n');
                fprintf(' * Instance  = %d ',i1);
                fprintf(' * Domain  = %d ',i2-1);
                fprintf(' * Iteration  = %d ',iteration);

                %seletion
                array = [fit(1) fit(2) fit(3) fit(4)];
                [sorted,index] = sort(array,'ascend');
             
                %crossover
                x(index(3)) = x(index(1));
                y(index(3)) = y(index(2));
                x(index(4)) = x(index(2));
                y(index(4)) = y(index(1));

                %mutation
                rand_num = randi([1 2],1,1);
                if rand_num == 1
                    x(index(4)) = (max_x - min_x)*rand(1,1) + min_x;
                else
                    y(index(4)) = (max_y - min_y)*rand(1,1) + min_y;
                end

                %fitness evaluation
                for i3 = 1:4
                    
                    Output_Parameters(i3,4)  = x(i3);
                    Output_Parameters(i3,6)  = y(i3);
                    
                    %eval ROM
                    [dofID,output] = CSMt_Identify_ROM(MESH, FE_SPACE, DATA, Output_Parameters(i3,:), node_positions, node_directions, ...
                    dofID, ROM_displ.V, Param_case(:,i1), figura, str_out);
            
                    fit(i3) = 0.5 * norm(observed.displ - output.displ,'fro')^2; %+ np.linalg.norm(x)
                
                end
                iteration = iteration + 1;
            end

            %brute force optimization for the subdomain
            if i2 == 1
                best_fit = fit(index(1));
                best_x = x(index(1));
                best_y = y(index(1));
                best_omega = i2-1;
            else
                if fit(index(1)) < best_fit
                    best_fit = fit(index(1));
                    best_x = x(index(1));
                    best_y = y(index(1));
                    best_omega = i2-1;
                end
            end    
        end
        dlmwrite([str_out,'/Frequency_identified.csv'],best_x,'-append')
        dlmwrite([str_out,'/Damage_level_identified.csv'],best_y,'-append')
        dlmwrite([str_out,'/Damage_class_identified.csv'],best_omega,'-append')
    end
end

