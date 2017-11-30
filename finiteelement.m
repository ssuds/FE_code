function [g_buckling_norm, stress_norm, dgbuckling_dAj, dsigma_dA] = finiteelement(design_variables,Nvars,file_name)
%% Read in input deck and store its values
[Nnodes,n_dof,nodes,Nmaterials,materials,Nelements,elements,NconcM,concmasses,Nforcecases,forcecases,NBCsets,BCdofID,Nstatic,static_load,Nmasscases,masscases,Ndynamic,dynamic_cases,sensitivityflag,~,~] = inputs(file_name);
 
%% Replace nominal design values with optimizer design values
if isempty(design_variables)
    display('Optimizer is being run with an empty input argument. Quantities will be computed for nominal design variable values.') 
else
    display('Optimizing against design variables specified in optimizer input deck.') 
    %Optimization has not been implemented yet.
end
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
stress = solvestress(Nelements,Nstatic,n_dof,c,forcecases,BCdofID,elements,nodes,dof_index,materials,u);

disp('The stresses in each element in PSI are displayed below in order of element number (vertically) and load case (horizontally)');
disp(stress);

%% Generate mass matrix for the unconstrained problem for each mass case 
Mglobal_masscases = MwithCases(Nmasscases,Mglobal,masscases,concmasses,dof_index);

%% Solve the Dynamic Cases to find Natural Frequencies, Eigenvectors, and Eigenvalues
[natural_freqs, phi, lambda] = solvedynamic(Ndynamic,dynamic_cases,Kglobal_bcs,Mglobal_masscases,n_dof);

%% Create buckling constraints for truss members
g_buckling = buckling(n_dof,Nelements,elements,nodes,dof_index,materials,Nstatic,stress);

%% Normalize the objective function and buckling constraints
stress_norm = normc(stress);
g_buckling_norm = normc(g_buckling);

%% Compute analytic sensitivities
if sensitivityflag == 1 | nargout >2 %Analytic sensitivities for will only be computed if the sensitivity flag is 1 within the input deck, or the optimizer code calling this function requests analytical sensitivities
    [dK_dA,dsigma_dA,dF_dA,dM_dA,dlambda_dA,dgbuckling_dAj] = sensitivities(n_dof,Nelements,elements,nodes,dof_index,materials,Kellocal_elements,T,Nstatic,u,static_load,Kglobal_bcs,c,Melglobal_elements,Ndynamic,phi,stress);
end

end