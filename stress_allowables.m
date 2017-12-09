function g_stress_allowables = stress_allowables(n_dof,Nelements,elements,nodes,dof_index,materials,Nstatic,stress)

[~,~,~,~,~,~,~,sigma_yt,sigma_yc] = elementproperties(n_dof,Nelements,elements,nodes,dof_index,materials);
g_stress_allowables = zeros(Nelements,Nstatic); %initialize matrix to hold buckling constraint values for all elements and static load cases

for ist = 1:Nstatic
    g_stress_allowables_case = zeros(Nelements,1); %initialize vector to hold buckling constraint values for each element
    for iel = 1:Nelements
        if stress(iel,ist) > 0 % Positive stresses are tensile
            g_stress_allowables_el = stress(iel,ist)-sigma_yt(iel,1); %If negative, the stress in the element is lower than the allowable. If positive, it is higher (the constraint is violated).
        elseif stress(iel,ist) < 0 % Negative stresses are compressive
            g_stress_allowables_el = -stress(iel,ist)-sigma_yc(iel,1); %If negative, the stress in the element is lower than the allowable. If positive, it is higher (the constraint is violated).
        elseif stress(iel,ist) == 0            
            g_stress_allowables_el = 0 - min(abs(sigma_yt(iel,1),abs(sigma_yc(iel,1)))); %This takes the conservative approach of choosing the lowest allowable stress for either compression or tension
        else
            disp('Check your stress allowable constraints');
        end
%           g_stress_allowables_el = abs(stress(iel,ist))-sigma_yt(iel,1); % Simplified version of the stress allowables constraint if the above fails. if the magnitude of stress exceeds the yield strength, this constraint will become positive (violated). Note - this assumes tensile and compressive yield strength are the same.
          g_stress_allowables_case(iel,1) = g_stress_allowables_el; %store the stress allowable constraint value for this element into the larger vector
    end
    g_stress_allowables(1:Nelements,ist) = g_stress_allowables_case; %store the stress allowable constraint values for this static load case into the larger matrix
end

end