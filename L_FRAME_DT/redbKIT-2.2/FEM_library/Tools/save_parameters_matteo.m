function save_parameters_matteo(str_out, param, index_param, index_out)

%     if index_param == 0
%         outfile = [str_out,'\param_used_sim_n_',num2str(index_out),str,'.txt'];
%     else
%         outfile = [str_out,'\param_used_sim_n_',num2str(index_param),str,'.txt'];
%     end
    outfile = [str_out,'\param_used_sim_n_',num2str(index_out),'_',num2str(index_param),'.txt'];
    if length(param) > 2
        i1 = 3;
    else
        i1 = 1;
    end
    datafile = fopen(outfile, 'wt');
    fprintf(datafile,'**************************************************************\n');
    fprintf(datafile,'Parameters employed\n');
    fprintf(datafile,'                - amplitude = %f\n',param(i1));
    fprintf(datafile,'                - frequency = %f\n',param(i1 + 1));
    
    fclose(datafile);
    
end