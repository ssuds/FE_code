%% Design Oriented Truss FE code for 2D structures
% AA 538 Introduction to Structural Optimization
% Shreyas Sudhakar

%% Clearing all previous data in workspace
clear; 

%% Define filename of input deck
file_name = 'input.txt';

%% Contents of input deck
type (file_name);

%% Read in input deck and store its values
[Nnodes,n_dof,nodes,Nmaterials,materials,Nelements,elements,NconcM,concmasses,Nforcecases,forcecases,NBCsets,BCdofID,Nstatic,static_load,Nmasscases,masscases,Ndynamic,dynamic_cases,sensitivityflag,Nvars,design_variables] = inputs(file_name);

i = 0;
while i < 10
    %% HW 3 Part 1
    [g_buckling, stress, dgbuckling_dAj, dsigma_dA] = finiteelement(design_variables, Nvars, file_name); %call the finite element code, provide the vector of design variables and recieve normalized vectors of constraints and the objective function, and also the analytic sensitivities of these quantities
    optimizer(g_buckling, stress, dgbuckling_dAj, dsigma_dA, file_name);
    i = i + 1;
end