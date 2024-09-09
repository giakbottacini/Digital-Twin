% Plot energy content POD

function plot_snapshot_energy(DATA,Sigma,S,V,which_type, i, pb_ID)

    %Sigma: autovalori ottenuti con la POD
    %S    : size of the snapshot matrix
    %V    : base ottenuta con la POD
    
    if nargin < 7
        pb_ID = ['/damaged_',DATA.pb_ID];
    end
    

    figure_output = figure;
    loglog(Sigma./Sigma(1),'-or');
    hold on
    cut_line_x = [1 : size(V,2)];
    cut_line_y = Sigma(size(V,2))./Sigma(1) * ones(1,size(V,2));
    loglog(cut_line_x, cut_line_y,'--k');
    grid on
    xlim([1  S])
    xlabel('\fontsize{16}number of snapshots')
    ylabel('\fontsize{16}error')
      
    if isempty(i)
        save_path = ['/Users/matteotorzoni/Desktop/Corigliano/Mensola_matteo/Dati/',DATA.problem,'/Danneggiamento/damaged_', DATA.pb_ID_f];
        saveas(figure_output,fullfile(save_path,[which_type]),'fig')
        saveas(figure_output,fullfile(save_path,[which_type]),'png')
    else
        [~, ~, ~] = mkdir(['/Users/matteotorzoni/Desktop/Corigliano/Mensola/Dati/Non_linear_fitting/Damage_Identification/',DATA.problem,'/damaged_',DATA.pb_ID_f, pb_ID], ['training_mu_',num2str(i)]);
        save_path = ['/Users/matteotorzoni/Desktop/Corigliano/Mensola/Dati/Non_linear_fitting/Damage_Identification/',DATA.problem,'/damaged_', DATA.pb_ID_f, pb_ID,'/training_mu_',num2str(i)];
        saveas(figure_output,fullfile(save_path,[which_type,'_training_param_',num2str(i)]),'fig')
        saveas(figure_output,fullfile(save_path,[which_type,'_training_param_',num2str(i)]),'png')
    end
                            
end                            

