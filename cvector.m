function c = cvector(n_dof,Nelements,elements,nodes,dof_index,materials)
% CVECTOR  Computes the c vectors for each element in a truss structure based on the number of degrees of freedom n_dof and the dof properties in dof_index.
%

c = zeros(Nelements,n_dof);%initialize empty c matrix
[~,theta,~,~,~,~,~] = elementproperties(n_dof,Nelements,elements,nodes,dof_index,materials); %Pull the angle of all elements in the structure

for iel=1:Nelements
    %find ids of node 1 and node 2 that the element connects from elements matrix
    node_1 = elements(iel,1);
    node_2 = elements(iel,2);

    %Identify the 2 dofs associated with node 1 of the element
    dof_node1_index = find(dof_index(:,2)==node_1); %find all instances of node_1 id within the node column of the dof_index matrix
    dof1_node1 = dof_index(dof_node1_index(1),1); %identify first dof associated with node 1
    dof2_node1 = dof_index(dof_node1_index(2),1); %identify second dof associated with node 1
    %Identify the 2 dofs associated with node 2 of the element
    dof_node2_index = find(dof_index(:,2)==node_2); %find all instances of node_2 id within the node column of the dof_index matrix
    dof1_node2 = dof_index(dof_node2_index(1),1); %identify first dof associated with node 2
    dof2_node2 = dof_index(dof_node2_index(2),1); %identify second dof associated with node 2

    %create temporary c vector
    c_vectemp = zeros(1,n_dof);
    %populate the c vector according to the element id to be able to calculate appropriate stresses for this element 
    %populate the first trig term in the c vector
    c_vectemp(1,dof1_node1)=-cos(theta(iel,1));
    %populate the second trig term in the c vector
    c_vectemp(1,dof2_node1)=-sin(theta(iel,1));
    %populate the third trig term in the c vector
    c_vectemp(1,dof1_node2)=cos(theta(iel,1));
    %populate the fourth trig term in the c vector
    c_vectemp(1,dof2_node2)=sin(theta(iel,1));
    %write temporary c vector into overall c matrix
    c(iel,1:n_dof) = c_vectemp;

end

end