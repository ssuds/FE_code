function [Kglobal_bcs,Mglobal_bcs] = KandMwithBCs(NBCsets,Kglobal,Mglobal,BCdofID,n_dof)
% Function creates global stiffness and mass matrices for each boundary conditions (BC) set
%%
%HW1 Step 5
Kglobal_bcs = cell(NBCsets,1); %initialize a cell array to store the global stiffness matrix for each set of boundary conditions
Kglobal_backup = Kglobal; %initialize a backup global stiffness matrix 
Mglobal_bcs = cell(NBCsets,1); %initialize a cell array to store the global mass matrix for each set of boundary conditions
Mglobal_backup = Mglobal; %initialize a backup global mass matrix 
for ibc=1:NBCsets
 Kglobal = Kglobal_backup; %Restore the backup of the stiffness matrix so that each boundary condition can be applied on a virgin stiffness matrix
 Mglobal = Mglobal_backup; %Restore the backup of the mass matrix so that each boundary condition can be applied on a virgin mass matrix 
 temp_boundary=BCdofID(1:n_dof,ibc); %store the vector of the current boundary condition set
 for i=1:n_dof
     diag_elem=Kglobal(i,i); %store the diagonal element i,i of the global stiffness matrix 
     if temp_boundary(i)==0 
         Kglobal(i,1:n_dof)=0; %zero out the ith row of the stiffness matrix
         Kglobal(1:n_dof,i)=0; %zero out the ith column of the stiffness matrix
         Kglobal(i,i)=diag_elem; %restore the i,i diagonal element to the stiffness matrix
         Mglobal(i,1:n_dof)=0; %zero out the ith row of the mass matrix
         Mglobal(1:n_dof,i)=0; %zero out the ith column of the mass matrix
         Mglobal(i,i)=0.0000000001; %enter a small value of mass on the diagonal, leading to uncoupled natural frequencies of the dummy 1dof oscillators corresponding to constrained dof that will be high and out of range of interest (simple trick to deal with zero motion dofs)
     end
 end
 Kglobal_bcs{ibc}=Kglobal; %Store the global stiffness matrix corresponding to the ibc'th boundary condition set into the stiffness matrix cell array
 Mglobal_bcs{ibc}=Mglobal; %Store the global mass matrix corresponding to the ibc'th boundary condition set into the mass matrix cell array
end
clear i; %clears out the id variable so that it doesnt mess with indexing later
clear ibc; %clears out the id variable so that it doesnt mess with indexing later