%% Design Oriented Truss FE code for 2D structures
% AA 538 Introduction to Structural Optimization
% Shreyas Sudhakar

%% Clearing all previous data in workspace
clear; 

%% Find number of elements from input deck
[~,~,~,~,~,Nelements,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = inputs('input.txt');

%% Design Variables
[design_variables, Nvars] = optinputs('input_optimizer.txt', Nelements); %Read in the design variables from the optimizer input deck

%% HW 3 Part 1
[g_buckling, stress, dgbuckling_dAj, dsigma_dA] = finiteelement(design_variables, Nvars); %call the finite element code, provide the vector of design variables and recieve normalized vectors of constraints and the objective function, and also the analytic sensitivities of these quantities


%% Normalized buckling constraints
display(g_buckling);
%% Sensitivity of the buckling constraints
display(dgbuckling_dAj);
%% Normalized objective function (stress in this instance)
display(stress);
%% Sensitivity of the objective function
display(cell2mat(dsigma_dA));
