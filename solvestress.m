function stress_matrix = solvestress(Nelements,Nstatic,n_dof,c,forcecases,BCdofID,elements,nodes,dof_index,materials,u)

[l,~,E,~,~,~,~] = elementproperties(n_dof,Nelements,elements,nodes,dof_index,materials); %compute commonly used element properties 

%find stresses

stress_matrix = zeros(Nelements,Nstatic); %initialize stress matrix to hold stresses for each static load case and each element
for ist=1:Nstatic
    stress_vec = zeros(1,Nelements); %initialize stress vector to capture stresses for each element
    for iel=1:Nelements
        stress = (E(iel,1)/l(iel,1))*c(iel,1:n_dof)*u(1:n_dof,ist); %compute stress {c}*{u}
        stress_vec(iel) = stress; %save the stress for the iel element into the iel'th element of the stress_vec
    end
    stress_matrix(1:Nelements,ist) = stress_vec;
end
end