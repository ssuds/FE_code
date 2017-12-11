function u = solvedisplacement(n_dof,Nstatic,static_load,force_inputs,Kglobal_bcs)

%HW1 Step 7: Solve for displacements and stresses in all static load cases
%find displacements
u = zeros(n_dof,Nstatic); %initialize displacement matrix to hold displacements for each static load case
for ist=1:Nstatic
    ibc = static_load(ist,2); %find the boundary condition set index ibc
    iforce = static_load(ist,1); %find the force set index iforce
    force = force_inputs(1:n_dof,iforce); %store the modified force vector associated to the index
    Kglobal_set = Kglobal_bcs{ibc}; %store the global stiffness matrix corresponding to BC set number ibc 
    u_case = mldivide(Kglobal_set,force); %Solve [K]{u} = {F} for {u}, lsqr function is used since K is sparse
    u(1:n_dof,ist)=u_case; %store the displacement for this load case in the Nstatic column of the displacement matrix
end

end