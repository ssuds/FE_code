function [Mglobal,Melglobal_elements] = massmatrix(n_dof,Nelements,elements,nodes,dof_index,materials)

%% Creation of global structural mass matrix

[l,theta,E,A,V,rho,m] = elementproperties(n_dof,Nelements,elements,nodes,dof_index,materials); %compute commonly used element properties 

%% Allocate memory for the global mass matrix for the structure: [Mglobal], initiate these matrices by setting all their elements to zero.

%initialize global mass matrix
Mglobal = zeros(n_dof,n_dof);
%section the global mass matrix into 2x2 sections
clear temp_ar;
temp_ar(1:(n_dof/2))=2;
Mglobal_sectioned = mat2cell(Mglobal,temp_ar,temp_ar);

%%
Melglobal_elements = cell(1,Nelements); %initialize cell matrix to store local mass matrix for each element
for iel=1:Nelements
    %find ids of node 1 and node 2 that the element connects from elements matrix
    node_1 = elements(iel,1);
    node_2 = elements(iel,2);
    
    
    %initialize an element mass matrix
    Melglobal = zeros(4,4);

    %adding the element mass to the diagonal of the element mass matrix,
    %corresponding to the 4 DOFs for n1 and n2 that the element connects
    Melglobal = diag([1 1 1 1]);
    Melglobal = m(iel,1)/2.*Melglobal;
    Melglobal_elements{iel} = Melglobal; %Store the element mass matrix for this element into the cell array of all elements
    %Partition [Melglobal] into 4 2 x 2 matrices
    Melglobal_sectioned = mat2cell(Melglobal,[2 2],[2 2]);
    %adding the element mass matrix into the global mass matrix in correct location
    Mglobal_sectioned{node_1,node_1} = Mglobal_sectioned{node_1,node_1} + Melglobal_sectioned{1,1};
    Mglobal_sectioned{node_1,node_2} = Mglobal_sectioned{node_1,node_2} + Melglobal_sectioned{1,2};
    Mglobal_sectioned{node_2,node_1} = Mglobal_sectioned{node_2,node_1} + Melglobal_sectioned{2,1};
    Mglobal_sectioned{node_2,node_2} = Mglobal_sectioned{node_2,node_2} + Melglobal_sectioned{2,2};


end
clear iel; %clears out the id variable so that it doesnt mess with indexing later


%Rebuild Mglobal from sectioned matrix
Mglobal = cell2mat(Mglobal_sectioned);
end