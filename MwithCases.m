function Mglobal_masscases = MwithCases(Nmasscases,Mglobal,masscases,concmasses,dof_index)


%% HW2 part a) Add mass matrix generation

Mglobal_masscases = cell(Nmasscases,1); %initialize a cell array to store the global mass matrix for each set of concentrated mass cases
Mglobal_backup = Mglobal; %initialize a backup global mass matrix 
for i=1:Nmasscases 
    for j=1:length(masscases{i})
        Mglobal = Mglobal_backup; %Restore the backup of the mass matrix so that each boundary condition can be applied on a virgin mass matrix 
        %find the mass value 
        Mj = concmasses(masscases{i}(j),1);
        %find the node the mass is located at
        nodej = concmasses(masscases{i}(j),2);
        %find the associated degrees of freedom for the node
        dof_node1_index = find(dof_index(:,2)==nodej); %find all instances of nodej id within the node column of the dof_index matrix
        i1 = dof_index(dof_node1_index(1),1); %identify first dof associated with node
        i2 = dof_index(dof_node1_index(2),1); %identify second dof associated with node
        %add the concentrated masses to the appropriate DOFs
        Mglobal(i1,i1)=Mglobal(i1,i1)+Mj;
        Mglobal(i2,i2)=Mglobal(i2,i2)+Mj;
    end 
 Mglobal_masscases{i}=Mglobal; %Store the global mass matrix corresponding to the i'th mass case set into the mass case matrix cell array
 end
clear i; %clears out the id variable so that it doesnt mess with indexing later
clear j; %clears out the id variable

end