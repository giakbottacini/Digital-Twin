function param_correct = CSM_correct_param(DATA,dim,param)


    %% Lamé constants computation
    % planestress/planestrain correction
    if dim == 2
        if contains(DATA.plane_analysis,'planestress')
            param(1,:)   = param(1,:)  .*  (1 + 2.*param(2,:)) ...
                                       ./ ((1 + param(2,:)).^2);
            param(2,:)   = param(2,:)  ./ (1 + param(2,:));
    
        elseif contains(DATA.plane_analysis,'planestrain')
            param(1,:)   = param(1,:) ;
            param(2,:)   = param(2,:);    
        end
    end

    % Lamé constants
        % 1 - mu coefficient
        % 2 - lambda coefficient
    param_correct      = param;
    param_correct(1,:) = param(1,:) / (2 + 2 .* param(2,:));
    param_correct(2,:) = param(1,:) .* param(2,:) ...
                                    ./ ((1 + param(2,:)).*(1 - 2 .* param(2,:)));


end