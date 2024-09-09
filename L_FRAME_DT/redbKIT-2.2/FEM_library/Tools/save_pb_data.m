function save_pb_data(DATA,mesh_name,str_sfx,M,N1,N2,dim,fem,node_position,...
                      direction_node,nodeID,compute_ei,vertices)
                  
    if isempty(str_sfx)
        current_folder=pwd;
        str_out=[current_folder, '\Graphs\'];
        outfile = [str_out,'data.txt'];
        datafile = fopen(outfile, 'wt');
        title='>>>>>>>>>>>  Elastic analysis data\n';
    else
        [~, ~, ~] = mkdir(['C:\Users\Rosafalco\Documents\Polimi\3_Dottorato\Dati\Non_linear_fitting\Damage_Identification\mensola\damaged_',...
                          DATA.pb_ID_f,'\'], ['damaged_', DATA.pb_ID]);
        str_out=['C:\Users\Rosafalco\Documents\Polimi\3_Dottorato\Dati\Non_linear_fitting\Damage_Identification\mensola\damaged_',DATA.pb_ID_f,'\damaged_', DATA.pb_ID];
        outfile = [str_out,'data.txt'];
        datafile = fopen(outfile, 'wt');
        title='>>>>>>>>>>>  Non linear fitting data_%s\n';
    end
    

    fprintf(datafile,'**************************************************************\n');
    fprintf(datafile,'Reduced Order Modelling of elasto-structural dynamics problems\n');
    fprintf(datafile,'\n');
    fprintf(datafile,title,str_sfx);
    fprintf(datafile,'\n');

    fprintf(datafile,'Problem dimension = %u\n',dim);
    if dim==2
        fprintf(datafile,'     ----> %s\n',DATA.plane_analysis);
    end
    fprintf(datafile,'FE type employed: "%s"\n',fem);
    
    fprintf(datafile,'energy content discrepancy of the POD model:\n');
    fprintf(datafile,'                - tol_POD_U = %.2e\n',DATA.tol_POD_U);
    fprintf(datafile,'parameters for the convergence of the reduced modes in time:\n');
    fprintf(datafile,'                - min_snapshots = %u --> minimum n° of snapshots required\n',DATA.min_snapshots);
    fprintf(datafile,'                - tol_P_k_k1_U  = %.2e --> |E(k+1)-E(k)|<P\n',DATA.tol_P_k_k1_U);
    fprintf(datafile,'                - tol_D_k_k1_U  = %.2e --> dev_std|E(k)|<D\n',DATA.tol_D_k_k1_U);
    fprintf(datafile,'          --> n° of selected modes = %u\n',M);
    fprintf(datafile,'                - tol_POD_Fext = %.2e\n',DATA.tol_POD_Fext);
    fprintf(datafile,'                - tol_POD_Fint = %.2e\n',DATA.tol_POD_Fint);
    fprintf(datafile,'          --> n° of selected points Fint = %u\n',N1);
    fprintf(datafile,'          --> n° of selected points Fint = %u\n',N2);
    fprintf(datafile,'\n');
    
    fprintf(datafile,'employed geometry and mesh discretization:\n');
    fprintf(datafile,'                - %s\n',mesh_name);
    fprintf(datafile,'\n');
    
    fprintf(datafile,'material features:\n');
    fprintf(datafile,'                - Material model = %s\n',DATA.Material_Model);
    if contains(DATA.Material_Model,'Linear') || contains(DATA.Material_Model,'StVenantKirchhoff')
        fprintf(datafile,'                - Young modulus = %.3e\n',DATA.Young);
        fprintf(datafile,'                - Poisson coefficient = %f\n',DATA.Poisson);
    end
    fprintf(datafile,'                - Density = %.3f\n',DATA.Density);
    fprintf(datafile,'\n');
    
    fprintf(datafile,'analysis time features:\n');
    fprintf(datafile,'                - time step = %.3e\n',DATA.time.dt);
    fprintf(datafile,'                - starting time = %.2f\n',DATA.time.t0);
    fprintf(datafile,'                - maximum training interval = %.2f\n',DATA.time.tf_short);
    fprintf(datafile,'                - analysis interval = %.2f\n',DATA.time.tf_long);
    fprintf(datafile,'\n');

    fprintf(datafile,'integration time features:\n');
    fprintf(datafile,'                - gamma = %.1f\n',DATA.time.gamma);
    fprintf(datafile,'                - beta = %.1f\n',DATA.time.beta);
    fprintf(datafile,'                - alpha_m = %.1f\n',DATA.time.alpha_m);
    fprintf(datafile,'                - alpha_f = %.1f\n',DATA.time.alpha_f);
    fprintf(datafile,'\n');
    
    fprintf(datafile,'introduced damage:\n');
    if DATA.damaged
        fprintf(datafile,'                - is damaged? --> yes\n');
        fprintf(datafile,'                - where?      --> subdomain n_%u',DATA.where_damaged);
    else
        fprintf(datafile,'                - is damaged? --> no\n');
    end

    if ~isempty(nodeID)
        fprintf(datafile,'\n');
        fprintf(datafile,'physical quantity inspected  --> displacement\n');
        for i1=1:length(nodeID)
            if direction_node(i1)==0
                direct='x';
            elseif direction_node(i1)==1
                direct='y';
            else
                direct='z';
            end
            if dim==2
                fprintf(datafile,'monitored node: %u in direction: %s (corresponding to coordinates x=%f y=%f)\n',...
                        nodeID(i1), direct, node_position(1,i1), node_position(2,i1));
            elseif dim==3
                fprintf(datafile,'monitored node: %u in direction: %s (corresponding to coordinates x=%f y=%f z=%f\n',...
                        nodeID(i1), direct, node_position(1,i1), node_position(2,i1),  node_position(3,i1));
            end
        end
    end
    
    if compute_ei
        [n, omega, f, T]=compute_eigen(DATA,vertices);
        fprintf(datafile,'\n');
        fprintf(datafile,'eigenfrequencies computed from an analytical formulation: %u\n',n);
        for i1=1:n
            fprintf(datafile,'     omega=%f \n', omega(i1));
            fprintf(datafile,'     f    =%f \n',f(i1));
            fprintf(datafile,'     T    =%T \n',T(i1));
        end
    end
    
    fclose(datafile);
    
end
