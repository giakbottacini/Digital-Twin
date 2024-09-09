% Plot energy content POD

function plot_snapshot_energy_thermo(DATA,Sigma,S,V,which_type, i)

    %Sigma: autovalori ottenuti con la POD
    %S    : size of the snapshot matrix
    %V    : base ottenuta con la POD
    

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
        [~, ~, ~] = mkdir(['/Users/matteotorzoni/Desktop/Corigliano/PORTALE/Dati/',DATA.problem,'/Temperatura/Temperature_',num2str(V)]);
        save_path = ['/Users/matteotorzoni/Desktop/Corigliano/PORTALE/Dati/',DATA.problem,'/Temperatura/Temperature_',num2str(V)];
        saveas(figure_output,fullfile(save_path,[which_type]),'fig')
        saveas(figure_output,fullfile(save_path,[which_type]),'fig')
    else
        [~, ~, ~] = mkdir(['/Users/matteotorzoni/Desktop/Corigliano/PORTALE/Dati/',DATA.problem,'/Temperatura/Temperature_',num2str(i)]);
        save_path = ['/Users/matteotorzoni/Desktop/Corigliano/PORTALE/Dati/',DATA.problem,'/Temperatura/Temperature_',num2str(i)];
        saveas(figure_output,fullfile(save_path,[which_type,'_training_param_',num2str(i)]),'png')
        saveas(figure_output,fullfile(save_path,[which_type,'_training_param_',num2str(i)]),'fig')
    end
                            
end                            

