%% Design Oriented Truss FE code for 2D structures
% AA 538 Introduction to Structural Optimization
% Shreyas Sudhakar

%% Clearing all previous data in workspace
clear; 

%% Define filename of input deck
file_name = 'input.txt';

%% Read in input deck and store its values
[Nnodes,n_dof,nodes,Nmaterials,materials,Nelements,elements,NconcM,concmasses,Nforcecases,forcecases,NBCsets,BCdofID,Nstatic,static_load,Nmasscases,masscases,Ndynamic,dynamic_cases,sensitivityflag] = inputs(file_name);

prompt = 'What mode would you like to run in?'; %0 to display the contents of the input deck, 1 for Finite Element Analysis of nominal values, 2 for constrained optimization to minimize mass in Brute Force mode
mode = input(prompt)

if mode == 1
    stress = finiteelement(file_name,[]);
elseif mode == 2
    A_opt = optimizer(file_name, [], []);
elseif mode == 0
    %% Contents of input deck
    type (file_name);
else
    disp('This is not a valid mode! Please run again with correct input.');
end