%% Truss FE code for 2D structures
% AA 538 Introduction to Structural Optimization
% Shreyas Sudhakar

%% Clearing all previous data in workspace
clear; 

%% Define filename of input deck
file_name = 'input.txt';

%% Contents of input deck
type (file_name);

%% Read in input deck and store its values
[Nnodes,n_dof,nodes,Nmaterials,materials,Nelements,elements,NconcM,concmasses,Nforcecases,forcecases,NBCsets,BCdofID,Nstatic,static_load,Nmasscases,masscases,Ndynamic,dynamic_cases,sensitivityflag] = inputs(file_name);
 
%% Defines the degrees of freedom and associates each dof with a node and a direction
dof_index = dofproperties(n_dof);

%% Compute the C vectors for each element
c = cvector(n_dof,Nelements,elements,nodes,dof_index,materials);

%% Compute the transformation matrix for each element
T = transformationmatrix(n_dof,Nelements,elements,nodes,dof_index,materials);

%% Create global stiffness matrix
[Kglobal,Kellocal_elements] = stiffnessmatrix(n_dof,Nelements,elements,nodes,dof_index,materials,T);

%% Create global mass matrix
[Mglobal,Melglobal_elements] = massmatrix(n_dof,Nelements,elements,nodes,dof_index,materials);

%% Create stiffness and mass matrices for each set of boundary conditions
[Kglobal_bcs,Mglobal_bcs] = KandMwithBCs(NBCsets,Kglobal,Mglobal,BCdofID,n_dof);

%% Create the force input vectors for each static load case
force_inputs = forceinputs(n_dof,Nstatic,static_load,forcecases,BCdofID);

%% Solve the static problem [K]{u}={F}
u = solvedisplacement(n_dof,Nstatic,static_load,force_inputs,Kglobal_bcs);

%% Compute and display the static solution (Stresses in each element)
stress_matrix = solvestress(Nelements,Nstatic,n_dof,c,forcecases,BCdofID,elements,nodes,dof_index,materials,u);

disp('The stresses in each element in PSI are displayed below in order of element number (vertically) and load case (horizontally)');
disp(stress_matrix);

%% Generate mass matrix for the unconstrained problem for each mass case 
Mglobal_masscases = MwithCases(Nmasscases,Mglobal,masscases,concmasses,dof_index);

%% Solve the Dynamic Cases to find Natural Frequencies, Eigenvectors, and Eigenvalues
[natural_freqs, phi, lambda] = solvedynamic(Ndynamic,dynamic_cases,Kglobal_bcs,Mglobal_masscases,n_dof);

%% Create buckling constraints for truss members
g_buckling = buckling(n_dof,Nelements,elements,nodes,dof_index,materials,Nstatic,stress_matrix);

%% Compute analytic sensitivities
if sensitivityflag == 1 %Analytic sensitivities for will only be computed if the sensitivity flag is 1 within the input deck
    [dK_dA,dsigma_dA,dF_dA,dM_dA,dlambda_dA,dgbuckling_dAj] = sensitivities(n_dof,Nelements,elements,nodes,dof_index,materials,Kellocal_elements,T,Nstatic,u,static_load,Kglobal_bcs,c,Melglobal_elements,Ndynamic,phi,stress_matrix);
end