function g_buckling = buckling(n_dof,Nelements,elements,nodes,dof_index,materials,Nstatic,stress_matrix)


%% HW2 Part e) Truss member buckling
[l,~,E,A,~,~,~] = elementproperties(n_dof,Nelements,elements,nodes,dof_index,materials);
g_buckling = zeros(Nelements,Nstatic); %initialize matrix to hold buckling constraint values for all elements and static load cases
for ist = 1:Nstatic
    g_buckling_case = zeros(Nelements,1); %initialize vector to hold buckling constraint values for each element
    for iel = 1:Nelements
      g_buckling_el = -1-stress_matrix(iel,ist)*((4*l(iel,1)^2)/(pi*E(iel,1)*A(iel,1))); %compute buckling criterion (if this is less than or equal to zero, buckling occurs)
      g_buckling_case(iel,1) = g_buckling_el; %store the buckling constraint value for this element into the larger vector
    end
    g_buckling(1:Nelements,ist) = g_buckling_case; %store the buckling constraint values for this static load case into the larger matrix
end

end