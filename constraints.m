function [g, g_eq] = constraints(file_name,A_opt,ibuckling)

%% Read in input deck and store its values
[Nnodes,n_dof,nodes,Nmaterials,materials,Nelements,elements,NconcM,concmasses,Nforcecases,forcecases,NBCsets,BCdofID,Nstatic,static_load,Nmasscases,masscases,Ndynamic,dynamic_cases,sensitivityflag] = inputs(file_name);

%% Defines the degrees of freedom and associates each dof with a node and a direction
dof_index = dofproperties(n_dof);

%% Call finite element code to compute stress 
stress = finiteelement(file_name,A_opt);

%% Build a matrix with all constraints
if ibuckling == 0 % Buckling constraints are not desired, return only stress constraints
    %% Create stress constraints
    g_stress_allowables = stress_allowables(n_dof,Nelements,elements,nodes,dof_index,materials,Nstatic,stress);
    %% Normalize the constraints
    g_stress_allowables_norm = normc(g_stress_allowables);
    %% Reshape constraints into a vector
    g = g_stress_allowables(:); %Reshape the constraints matrix into a vector
elseif ibuckling == 1 % Buckling constraints are desired, return both stress and buckling constraints in one constraint vector
    %% Create stress constraints
    g_stress_allowables = stress_allowables(n_dof,Nelements,elements,nodes,dof_index,materials,Nstatic,stress);
    %% Create buckling constraints for truss members
    g_buckling = buckling(n_dof,Nelements,elements,nodes,dof_index,materials,Nstatic,stress);
    %% Normalize the constraints
    g_stress_allowables_norm = normc(g_stress_allowables);
    g_buckling_norm = normc(g_buckling);
    %% Assembly of g vector
    g = horzcat(g_stress_allowables_norm,g_buckling_norm);
    g = g(:) %Reshape the constraints matrix into a vector
else
    disp('Undefined input to constraints solver, check the constraints');
    g = [];
end
g_eq = []; %We have no equality constraints!

end
