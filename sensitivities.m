 function [dK_dA,dsigma_dA,dF_dA,dM_dA,dlambda_dA,dgbuckling_dAj] = sensitivities(n_dof,Nelements,elements,nodes,dof_index,materials,Kellocal_elements,T,Nstatic,u,static_load,Kglobal_bcs,c,Melglobal_elements,Ndynamic,phi,stress)

 [l,theta,E,A,V,rho,m] = elementproperties(n_dof,Nelements,elements,nodes,dof_index,materials); %Compute commonly used properties for all elements 
%% HW2 part c) Analytic sensitivities with respect to cross sectional areas - static cases

%initialize dK_dA cell matrix, which will hold the d[K]/dA used in computing analytic sensitivities for each element
dK_dA = cell(Nelements,1);

clear temp_ar;
temp_ar(1:(n_dof/2))=2; 
for iel = 1:Nelements
    %find ids of node 1 and node 2 that the element connects from elements matrix
    node_1 = elements(iel,1);
    node_2 = elements(iel,2);
    dK_dAglobal = zeros(n_dof,n_dof); %initialize the global d[K]/dA matrix which is used for analytic sensitivities
    dK_dAglobal_sectioned = mat2cell(dK_dAglobal,temp_ar,temp_ar); %Section dK/dA matrix into 2x2 sections
    %build the element contribution to the global dK/dA matrix
    dK_dAlocal = zeros(4,4);
    dK_dAlocal(1,1) = 1;
    dK_dAlocal(3,1) = -1;
    dK_dAlocal(1,3) = -1;
    dK_dAlocal(3,3) = 1;
    dK_dAlocal = (E(iel,1)/l(iel,1)).*Kellocal_elements{iel};
    dK_dAelglobal = transpose(T{iel})*dK_dAlocal*T{iel}; %transform the element dK/dA matrix to the global coordinate system
    dK_dAelglobal_sectioned = mat2cell(dK_dAelglobal,[2 2],[2 2]); %partition the dK/dA for the element in the same way
    %position dK_dA element into global coordinates for the element DOFs
    dK_dAglobal_sectioned{node_1,node_1} = dK_dAelglobal_sectioned{1,1};
    dK_dAglobal_sectioned{node_1,node_2} = dK_dAelglobal_sectioned{1,2};
    dK_dAglobal_sectioned{node_2,node_1} = dK_dAelglobal_sectioned{2,1};
    dK_dAglobal_sectioned{node_2,node_2} = dK_dAelglobal_sectioned{2,2};

    %Rebuild dK_dAglobal from sectioned matrix
    dK_dAglobal = cell2mat(dK_dAglobal_sectioned);

    %Store the dK_dA in the overall dK_dA cell array corresponding to the element it is associated to 
    dK_dA{iel} = dK_dAglobal;

end

dsigma_dA = cell(Nstatic,1); %initialize stress sensitivity matrix to hold stress sensitivity for each static load case and each element
for ist=1:Nstatic
    %initialize variables to compute and store stress sensitivity for each element
    dF_dA = zeros(n_dof,1); %compute d{F}/dA - currently initialized as a zero vector as the forces being applied have no dependence on our structural design variables
    pseudo_load_loadcase = zeros(n_dof,Nelements); %initialize pseudo_load cell array to store pseudo load for each element under this load case
    du_dA_loadcase = zeros(n_dof,Nelements); %initialize d{u}/dA cell array to store d{u}/dA for each element under this load case
    dsigma_dA_loadcase = zeros(n_dof,Nelements); %initialize dsigma/dA cell array to store dsigma/dA for each element under this load case
    %compute stress sensitivity for element
    u_case = u(1:n_dof,Nstatic) %recall the computed displacements for this load case
    ibc = static_load(ist,2); %find the boundary condition set used in this static load case
    Kglobal_set = Kglobal_bcs{ibc}; %store the global stiffness matrix corresponding to BC set number ibc 
    for iel=1:Nelements
        c_element = c(iel,1:n_dof); %get the vector {c} for the element
        pseudo_load_loadcase(1:n_dof,iel) = dF_dA-(dK_dA{iel}*u_case); %compute pseudo load for each element
        du_dA_loadcase(1:n_dof,iel) = lsqr(Kglobal_set,pseudo_load_loadcase(1:n_dof,iel)); %Solve [dK/dA]{du/dA} = {pseudo_load} for {du/dA}, lsqr function is used since K is sparse
        dsigma_dA_loadcase(1:n_dof,iel) = c_element*du_dA_loadcase(1:n_dof,iel); %compute stress sensitivity for the element NOTE: this should be {c-transpose}*d{u}/dA but it appears one of these is a row vector when it should be a column vector
    end
    dsigma_dA{ist} = dsigma_dA_loadcase;
end
%% HW2 part d) Analytic sensitivities with respect to cross sectional areas - dynamic cases
%initialize dM_dA cell matrix, which will hold the d[M]/dA used in computing analytic sensitivities for each element
dM_dA = cell(Nelements,1);

clear temp_ar;
temp_ar(1:(n_dof/2))=2;  
for iel = 1:Nelements
    dM_dAglobal = zeros(n_dof,n_dof); %initialize the global d[M]/dA matrix which is used for analytic sensitivities
    dM_dAglobal_sectioned = mat2cell(dM_dAglobal,temp_ar,temp_ar); %Section dM/dA matrix into 2x2 sections 
    dM_dAelglobal = ((rho(iel,1)*l(iel,1))/2).*(Melglobal_elements{iel}); %Generate the element contribution to the global dM/dA matrix
    dM_dAelglobal_sectioned = mat2cell(dM_dAelglobal, [2 2],[2 2]);

    %adding the element mass sensitivity matrix into the global mass sensitivity matrix in correct location
    dM_dAglobal_sectioned{node_1,node_1} = dM_dAelglobal_sectioned{1,1};
    dM_dAglobal_sectioned{node_1,node_2} = dM_dAelglobal_sectioned{1,2};
    dM_dAglobal_sectioned{node_2,node_1} = dM_dAelglobal_sectioned{2,1};
    dM_dAglobal_sectioned{node_2,node_2} = dM_dAelglobal_sectioned{2,2};

    %Rebuild dMdAglobal from sectioned matrix
    dM_dAglobal = cell2mat(dM_dAglobal_sectioned);

    %Store the dK_dA in the overall dK_dA cell array corresponding to the element it is associated to 
    dM_dA{iel} = dM_dAglobal;
end


dlambda_dA = cell(1,Ndynamic); %Initialize a cell array to hold the dlambda/dA's computed for each dynamic case
for idyn = 1:Ndynamic
    dlambda_dA_case = cell(1,Nelements); %Initialize a cell array to hold the dlambda/dA computed for each element in this dynamic case
    phi_case = phi{idyn}; %retrieve the eigenvector for this dynamic case 
    for iel = 1:Nelements
        numerator = (phi_case.')*(dK_dA{iel}-(lambda_case*(dM_dA{iel})))*phi_case;
        denominator = (phi_case.')*Mglobal_set*phi_case;
        dlambda_dA_el = numerator/denominator; %calculate the analytic sensitivity for this element
        dlambda_dA_case{iel} = dlambda_dA_el; %store the computed analytic sensitivity 
    end

     dlambda_dA{idyn} = dlambda_dA_case; %store the analytic sensitivity computed for this case in the dlambda/dA cell array
end
clear idyn; %clears out the id variable so that it doesnt mess with indexing later
clear iel; %clears out the id variable so that it doesnt mess with indexing later

%% HW2 Part e) compute sensitivity of buckling constraints 
dgbuckling_dAj = zeros(Nelements); %initialize array to store sensitivities of buckling constraints in all elements with respect to areas of all elements            
for ist = 1:Nstatic
    for ieli = 1:Nelements
        for ielj = 1:Nelements
            if ielj == ieli %the second term in the sensitivity equation appears only when calculating the sensitivity of a buckling constraint in an element wrt the area of the same element
                dgbucklingi_dAj_ielj = (dsigma_dA{ist}(1,ielj)*((4*l(ieli,1)^2)/(pi*E(ieli,1)*A(ieli,1))))+(stress(ieli,ist)*((4*l(ieli,1)^2)/(pi*E(ieli,1)*(A(ielj,1)^2))));
            else
                dgbucklingi_dAj_ielj = dsigma_dA{ist}(1,ielj)*((4*l(ieli,1)^2)/(pi*E(ieli,1)*A(ieli,1)));
            end
            dgbuckling_dAj_ieli(ielj) = dgbucklingi_dAj_ielj; %store buckling constraints in one element with respect to areas of all elements
        end 
        dgbuckling_dAj(1:Nelements,ieli) = dgbuckling_dAj_ieli; %store buckling constraints in all elements with respect to areas of all elements
    end
end       
 end