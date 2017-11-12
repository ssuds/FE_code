function [Kglobal,Kellocal_elements] = stiffnessmatrix(n_dof,Nelements,elements,nodes,dof_index,materials,T)

%% Creation of global stiffness & structural mass matrices

[l,theta,E,A,V,rho,m] = elementproperties(n_dof,Nelements,elements,nodes,dof_index,materials); %compute commonly used element properties 

%% Allocate memory for the global stiffness matrix for the structure: [Kglobal] initiate these matrices by setting all their elements to zero.

%initialize global stiffness matrix
Kglobal = zeros(n_dof,n_dof);
%section the global stiffness matrix into 2x2 sections
clear temp_ar;
temp_ar(1:(n_dof/2))=2;
Kglobal_sectioned = mat2cell(Kglobal,temp_ar,temp_ar);

%%
% Step 3

Kellocal_elements = cell(1,Nelements); %initialize cell matrix to store local stiffness matrix for each element
for iel=1:Nelements
    %find ids of node 1 and node 2 that the element connects from elements matrix
    node_1 = elements(iel,1);
    node_2 = elements(iel,2);

    %generate the local stiffness matrix for the element
    %note that A is multiplied in after the fact to streamline the code when optimization is in play

    Kellocal = zeros(4,4);
    Kellocal(1,1) = 1;
    Kellocal(3,1) = -1;
    Kellocal(1,3) = -1;
    Kellocal(3,3) = 1;
    Kellocal = (E(iel,1)/l(iel,1)).*Kellocal;
    Kellocal = A(iel,1).*Kellocal;
    Kellocal_elements{iel} = Kellocal; %store the local stiffness matrix in a cell array

    %Find the 4 x 4 element stiffness matrix in global coordinates [KelGlobal]=[Ttranspose][KelLocal][T]
    Kelglobal = transpose(T{iel})*Kellocal*T{iel};
    %Partition [KelGlobal] into 4 2 x 2 matrices
    Kelglobal_sectioned = mat2cell(Kelglobal,[2 2],[2 2]);
    %Add element global stiffness matrix to global stiffness matrix in the
    %correct locations
    Kglobal_sectioned{node_1,node_1} = Kglobal_sectioned{node_1,node_1} + Kelglobal_sectioned{1,1};
    Kglobal_sectioned{node_1,node_2} = Kglobal_sectioned{node_1,node_2} + Kelglobal_sectioned{1,2};
    Kglobal_sectioned{node_2,node_1} = Kglobal_sectioned{node_2,node_1} + Kelglobal_sectioned{2,1};
    Kglobal_sectioned{node_2,node_2} = Kglobal_sectioned{node_2,node_2} + Kelglobal_sectioned{2,2};

end
clear iel; %clears out the id variable so that it doesnt mess with indexing later

%Rebuild Kglobal from sectioned matrix
Kglobal = cell2mat(Kglobal_sectioned);

end