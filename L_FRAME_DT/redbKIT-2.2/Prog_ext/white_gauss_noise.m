function [monitored_gdls] = white_gauss_noise(monitored_gdls, dev_std)

    %description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % white_gauss_noise.m
    %adds to the input signal a white gaussian noise (simple randomly
    %sampled)
    % Author: Luca Rosafalco
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    invcumfunct  = @(x)(norminv(x, 0, dev_std));                           % function handle - cumulative density function per probabilit? gaussiana
    
    for i1 = 1 : length(monitored_gdls)
        monitored_gdls(i1) = monitored_gdls(i1) + invcumfunct(rand(1));    % sampling compatibile con questa pdf
    end
    
end