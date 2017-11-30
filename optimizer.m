function [] = optimizer(g_buckling, stress, dgbuckling_dAj, dsigma_dA, file_name)
%% Find number of elements from input deck
[~,~,~,~,~,Nelements,~,~,~,~,~,~,~,~,~,~,~,~,~,~,design_variables, Nvars] = inputs(file_name);
%% Normalized buckling constraints
display(g_buckling);
%% Sensitivity of the buckling constraints
display(dgbuckling_dAj);
%% Normalized objective function (stress in this instance)
display(stress);
%% Sensitivity of the objective function
display(cell2mat(dsigma_dA));
end