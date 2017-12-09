function A_opt = optimizer(file_name, dgbuckling_dAj, dsigma_dA)
%% Read in input deck
[Nnodes,n_dof,nodes,Nmaterials,materials,Nelements,elements,NconcM,concmasses,Nforcecases,forcecases,NBCsets,BCdofID,Nstatic,static_load,Nmasscases,masscases,Ndynamic,dynamic_cases,sensitivityflag] = inputs(file_name);
[~,~,~,A,~,~,~] = elementproperties(n_dof,Nelements,elements,nodes,[],materials); % Find initial areas from input deck

%% Optimization with fmincon

min = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]; %Minimum allowable gauge for each element in the problem

[A_optimized,stress_optimized] = fmincon(@(A_opt)structure_mass(file_name,A_opt),A,[],[],[],[],min,[],@(A_opt)constraints(file_name,A_opt,0));

A_opt = A_optimized
end