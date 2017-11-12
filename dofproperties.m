function dof_index = dofproperties(n_dof)
% DOFPROPERTIES Defines the degrees of freedom and associates each dof with a node and a direction
% 

%initialize matrix dof_index that correlates each DOF to the corresponding node and direction
dof_index = zeros(n_dof,3);
%populate matrix dof_index
for i=1:n_dof
 dof_index(i,1)=i; %dof number
 dof_index(i,2)=ceil(i/2); %node number
 dof_index(i,3)=2-(mod(i,2)); %direction (1 for x, 2 for y)
end
clear i; %clears out the id variable so that it doesnt mess with indexing later

end