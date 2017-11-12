function T = transformationmatrix(n_dof,Nelements,elements,nodes,dof_index,materials)

[~,theta,~,~,~,~,~] = elementproperties(n_dof,Nelements,elements,nodes,dof_index,materials); %compute angle for each element 

T = cell(1,Nelements); %initialize cell matrix to store transpose matrix for each element
for iel=1:Nelements
%generate the transformation matrix for the element
T_element = zeros(4,4);
T_2x2 = zeros(2,2);
T_2x2(1,1)=(cos(theta(iel,1)));
T_2x2(1,2)=(sin(theta(iel,1)));
T_2x2(2,1)=(-sin(theta(iel,1)));
T_2x2(2,2)=(cos(theta(iel,1)));
% Section T into 2x2 cells so the T_2x2 matrices can be merged in
T_sectioned = mat2cell(T_element,[2 2],[2 2]);
%build the 4x4 T matrix by putting the T_2x2 matrices on the diagonal
T_sectioned{1,1} = T_2x2;
T_sectioned{2,2} = T_2x2;
%rebuild the T matrix from the sectioned T matrix
T_element=cell2mat(T_sectioned);
%store T matrix for each element
T{iel} = T_element;
end

end