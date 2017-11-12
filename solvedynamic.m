function [natural_freqs, phi, lambda] = solvedynamic(Ndynamic,dynamic_cases,Kglobal_bcs,Mglobal_masscases,n_dof)


%% HW2 part b) Define the number of natural frequency / mode shape cases
%% Dynamic Cases

natural_freqs = cell(1,Ndynamic); %initialize a cell array to store the natural frequencies found for each dynamic case
phi = cell(1,Ndynamic); %initialize a cell array to store the eigenvectors phi found for each dynamic case
lambda = cell(1,Ndynamic); %initialize a cell array to store the eigenvalues lambda found for each dynamic case

for idyn=1:Ndynamic
 ibc = dynamic_cases(idyn,1); %find the boundary condition set index ibc
 imcase = dynamic_cases(idyn,2); %find the mass case set id
 nlow = dynamic_cases(idyn,3); %find the number of lowest natural frequencies we want to concider, nlow
 Kglobal_set = Kglobal_bcs{ibc}; %store the global stiffness matrix corresponding to BC set number ibc (SPCs already applied)
 Mglobal_set = Mglobal_masscases{imcase}; %store the global mass matrix corresponding to mass set number imcase (no need to apply SPCs since matrix is diagonal and thus decouples SPC'ed dofs from the rest)

 % Hack to eliminate decoupled frequencies corresponding to SPCed dofs that are in the range of frequencies of the constrained structure (they do not represent frequencies of the constrained structure and are just a consequence of the trick by which we imposed SPCs.
 temp_boundary=BCdofID(1:n_dof,ibc); %store the vector of the current boundary condition set
 for i=1:n_dof
     diag_elem=Kglobal_set(i,i); %store the diagonal element i,i of the global stiffness matrix for this boundary condition set
     if temp_boundary(i)==0 
         Kglobal_set(i,i)=10000*Kglobal_set(i,i); %multiply the i,i diagonal element of the stiffness matrix by 10000 to move it out of our frequency range of interest
     end
 end
 
 [phi_case,lambda_case] = eig(Kglobal_set,Mglobal_set); %Solve the natural vibrations generalized eigenvalue problem lambda[M]{phi}=[K]{phi}
 lambdas_sorted = sort(lambda_case(lambda_case>0)); %Sort the lambda array from lowest to highest and eliminate zeros
 natural_freqs_sorted = (lambdas_sorted.^(1/2)).*(1/2*pi); %Convert the lambdas to natural frequencies in hertz
 natural_freqs_case = natural_freqs_sorted(1:nlow); %Choose the nlow lowest natural frequencies and store in a vector
 natural_freqs{idyn} = natural_freqs_case; %store the natural frequencies found for this case in the natural frequencies cell array
 lambda{idyn} = lambda_case; %store the eigenvalue found for this case in the eigenvalue cell array
 phi{idyn} = phi_case; %store the eigenvector found for this case in the eigenvector cell array
end
clear idyn; %clears out the id variable so that it doesnt mess with indexing later

end