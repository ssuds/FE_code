function force_inputs = forceinputs(n_dof,Nstatic,static_load,forcecases,BCdofID);


%HW1 Step 6: create the force input vectors for each static load case

force_inputs = zeros(n_dof,Nstatic); %initialize new matrix of force inputs for load cases 
force_new = zeros(n_dof,1); %initialize new force vector
for ist=1:Nstatic
    iforce = static_load(ist,1); %find the force vector index iforce
    ibc = static_load(ist,2); %find the boundary condition set index ibc
    force = forcecases(1:n_dof,iforce); %store the force vector associated to the index
    bc = BCdofID(1:n_dof,ibc); %store the boundary condition vector associated to the index
    for i=1:n_dof
        if bc(i) == 0 %check whether element in force vector corresponds to a dof that is constrained to zero motion
            force_new(i) = 0; % zero out the element in the new force vector
        else
            force_new(i) = force(i); % retain the element in the new force vector if there isn't a constraint
        end
    end
    force_inputs(1:n_dof,ist) = force_new; % create a new matrix of force inputs for load cases, where each column contains the forces applied to the loads case modified to account for the boundary conditions of the load case
end

end