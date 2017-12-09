function mass = structure_mass(file_name,A_opt)
%% Read in input deck 
[~,n_dof,nodes,~,materials,Nelements,elements,~,~,~,~,~,~,~,~,~,~,~,~,~] = inputs(file_name);
%% Defines the degrees of freedom and associates each dof with a node and a direction
dof_index = dofproperties(n_dof);
%% store vectors containing length and density of each element
[l,~,~,~,~,rho,~,~,~] = elementproperties(n_dof,Nelements,elements,nodes,dof_index,materials);
A_opt
mass_el = zeros(Nelements,1); % Initialize vector to hold element masses 
for iel = 1:Nelements
    mass_el(iel) = l(iel)*rho(iel)*A_opt(iel);
end

mass = sum(mass_el)
