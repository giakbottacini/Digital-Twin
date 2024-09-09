function [str_sfx]=save_for_ANN(U,U_POD,t)

    % general info
    str_out='C:\Users\Rosafalco\Documents\Polimi\3_Dottorato\Dati\Non_linear_fitting';
    str_prx='E_deter_const_';
    fprintf('Warning: the file will be saved with: \n')
    fprintf('              %s* \n',str_prx)
    fprintf('as initial string and with: \n')
    fprintf('%s \n',str_out)
    fprintf('as saving path \n')
    
    % require the user to specify the suffix of the saved data
    prompt='What is the suffix you want to use to save the data? e.g. 2 >>';
    str_sfx=num2str(input(prompt));
    
    % save the data of the FEM simulation
    U_output=zeros(size(U,1)+1,length(t));
    U_output(1,:)=t;
    U_output(2:size(U_output,1),:)=U(:,1:length(t));
    save([str_out,'\',str_prx,'U_',str_sfx,'.mat'],'U_output')
    dlmwrite([str_out,'\',str_prx,'U_',str_sfx,'.csv'],U_output,'delimiter',' ')
    
    % save the data of the ROM simulation
    U_POD_output=zeros(size(U_POD,1)+1,length(t));
    U_POD_output(1,:)=t;
    U_POD_output(2:size(U_POD_output,1),:)=U_POD(:,1:length(t));
    save([str_out,'\',str_prx,'U_POD_',str_sfx],'U_POD_output')
    dlmwrite([str_out,'\',str_prx,'U_POD_',str_sfx,'.csv'],U_POD_output,'delimiter',' ')
    
    % save the value pf the discrepancy between them
    U_discr=U-U_POD;
    U_discr_output=zeros(size(U_discr,1)+1,length(t));
    U_discr_output(1,:)=t;
    U_discr_output(2:size(U_discr_output,1),:)=U_discr(:,1:length(t));
    save([str_out,'\',str_prx,'U_discr_',str_sfx],'U_discr_output')
    dlmwrite([str_out,'\',str_prx,'U_discr_',str_sfx,'.csv'],U_discr_output,'delimiter',' ')
    
    
end