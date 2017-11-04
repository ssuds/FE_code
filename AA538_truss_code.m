%% Truss FE code for 2D structures - AA 538 Introduction to Structural Optimization - Shreyas Sudhakar

%Clearing all previous data in workspace
clear; 

%% Input file (This is stored externally as input.txt but shown here for reference)

%%
% 
% 
% 
% 
% 
%  #Number of nodes
%  6
%  #Coordinates of nodes
%  #id	x	y
%  1 720 360
%  2 720 0
%  3 360 360
%  4 360 0
%  5 0 360
%  6 0 0
%  #Number of materials
%  2
%  #Material Properties
% #id	E	sigma_yt	sigma_yc	rho
% 1 10000000 25000 25000 0.1
% 2 10000000 75000 75000 0.1
% #Number of elements
% 10
% #element properties
% #element_id	id_node1	id_node2	id_material	A
% 1 5 3 1 7.90
% 2 3 1 1 0.10
% 3 6 4 1 8.10
% 4 4 2 1 3.90
% 5 3 4 1 0.10
% 6 1 2 1 0.10
% 7 5 4 1 5.80
% 8 6 3 1 5.51
% 9 3 2 2 3.68
% 10 4 1 1 0.14
% #Number of concentrated non structural masses
% 1
% #concentrated mass properties
% #id	mass	node
% 1 0 1
% #Number of static force cases
% 1
% #force case properties
% #id	forces(zeros for unloaded dofs or non zero for forces in that direction and node)
% 1 0 0 0 -100000 0 0 0 -100000 0 0 0 0
% #Number of Boundary condition sets
% 1
% #boundary condition set properties
% #id	boundary_condition(free 1 or constrained 0)
% 1 1 1 1 1 1 1 1 1 0 0 0 0
% #Number of static load cases
% 1
% #static load case properties
% #id	forcevector_index	boundary_condition_case		
% 1 1 1
% #Number of natural vibrations cases
% 1
% #dynamic load case properties
% #id	boundary_condition_case	number_lowest_frequencies
% 1 1 2



%% Inputs for global stiffness and structural mass matrices

%Read input file, this is basically my punchcard!
fid=fopen('input.txt','r'); 
temp_str=fgetl(fid); %read over the first line, which is not used for anything other than formatting

%read in number of nodes as an integer
Nnodes=fscanf(fid,'%i',1); 

%read over the next three lines, the first two which are not used for
%anything other than formatting and third last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%read in the coordinates for each node and save in the array nodes
nodes=zeros(Nnodes,2); %initialize the nodes array, which has two columns for x and y
for i=1:Nnodes
 temp_ar=fscanf(fid,'%i %f %f',3); %reads in each node into a temporary array consisting of id, x coordinate, y coordinate
 nodes(temp_ar(1),1:2)=temp_ar(2:3); %assigns the x and y coordinates from the temporary array to the id'th row of the nodes matrix
end
clear i; %clears out the id variable so that it doesnt mess with indexing later

%read over the next two lines, the first which is not used for
%anything other than formatting and last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%read in number of materials as an integer
Nmaterials = fscanf(fid,'%i',1);

%read over the next three lines, the first two which are not used for
%anything other than formatting and third last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);
temp_str=fgetl(fid);

materials=zeros(Nmaterials,4); %initialize the materials array, which has 4 properties which we are tracking
for j=1:Nmaterials
 temp_ar=fscanf(fid,'%i %f %f %f %f',5); %reads in each material into a temporary array consisting of id, Young’s Modulus E(j), Yield stress in tension sigma_yt(j), Yield stress in compression sigma_yc(j), density rho(j), 
 materials(temp_ar(1),1:4)=temp_ar(2:5); %assigns the properties from the temporary array to the id'th row of the materials matrix
end
clear j; %clears out the id variable so that it doesnt mess with indexing later

%read over the next two lines, the first which is not used for
%anything other than formatting and last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%Read in Number of elements as an integer
Nelements = fscanf(fid,'%i',1);

%read over the next three lines, the first two which are not used for
%anything other than formatting and third last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);
temp_str=fgetl(fid);

elements=zeros(Nelements,4); %initialize the elements array, which has 4 properties which we are tracking
for k=1:Nelements
 temp_ar=fscanf(fid,'%i %f %f %f %f',5); %reads in each element into a temporary array consisting of element_id k from 1 to Nelements, the id numbers of the two nodes it connects: n_1(k),n_2(k), the material id for this element id_material(k), and the cross sectional area of the element A(k)
 elements(temp_ar(1),1:4)=temp_ar(2:5); %assigns the properties from the temporary array to the id'th row of the elements matrix
end
clear k; %clears out the id variable so that it doesnt mess with indexing later


%% Creation of global stiffness & structural mass matrices
% Step 1

%Define degrees of freedom (2 times the number of nodes for our 2d anal)
n_dof = Nnodes*2;
%initialize matrix dof_index that correlates each DOF to the corresponding 
%node and direction
dof_index = zeros(n_dof,3);
%populate matrix dof_index
for i=1:n_dof
 dof_index(i,1)=i; %dof number
 dof_index(i,2)=ceil(i/2); %node number
 dof_index(i,3)=2-(mod(i,2)); %direction (1 for x, 2 for y)
end
clear i; %clears out the id variable so that it doesnt mess with indexing later

%%
% Step 2

%initialize global stiffness matrix
Kglobal = zeros(n_dof,n_dof);
%section the global stiffness matrix into 2x2 sections
clear temp_ar;
temp_ar(1:(n_dof/2))=2;
Kglobal_sectioned = mat2cell(Kglobal,temp_ar,temp_ar);
%initialize global mass matrix
Mglobal = zeros(n_dof,n_dof);
%section the global mass matrix into 2x2 sections
clear temp_ar;
temp_ar(1:(n_dof/2))=2;
Mglobal_sectioned = mat2cell(Mglobal,temp_ar,temp_ar);
%%
% Step 3
%initialize empty c matrix
c = zeros(Nelements,n_dof);
for iel=1:Nelements
%find ids of node 1 and node 2 that the element connects from elements
%matrix
node_1 = elements(iel,1);
node_2 = elements(iel,2);
%find x and y coordinates of node 1 and node 2 from the nodes matrix
node_1x = nodes(node_1,1);
node_1y = nodes(node_1,2);
node_2x = nodes(node_2,1);
node_2y = nodes(node_2,2);
%calculate length of element 
l = sqrt((node_2y-node_1y)^2+(node_2x-node_1x)^2);
%calculate the angle the line from node 1 to node 2 creates with the
%global x axis
theta=atan2(node_2y-node_1y,node_2x-node_1x);
%Find the id_material for the element
id_material = elements(iel,3);
%Find material properties based on the id_material
E=materials(id_material,1);
rho=materials(id_material,4);
%Find the cross sectional area of the element A-iel
A = elements(iel,4);
%Find the volume of the element A-iel * L-iel (note: the –iel is an index
%not a subtraction)
V = A*l;
%Find the mass of the element m-iel=elementMaterialDensity*elementVolume
m=rho*V;
%generate the local stiffness matrix for the element
%note that A is multiplied in after the fact to streamline the code when
%optimization is in play
Kellocal = zeros(4,4);
Kellocal(1,1) = 1;
Kellocal(3,1) = -1;
Kellocal(1,3) = -1;
Kellocal(3,3) = 1;
Kellocal = (E/l).*Kellocal;
Kellocal = A.*Kellocal;
%generate the transformation matrix for the element
T = zeros(4,4);
T_2x2 = zeros(2,2);
T_2x2(1,1)=(cos(theta));
T_2x2(1,2)=(sin(theta));
T_2x2(2,1)=(-sin(theta));
T_2x2(2,2)=(cos(theta));
% Section T into 2x2 cells so the T_2x2 matrices can be merged in
T_sectioned = mat2cell(T,[2 2],[2 2]);
%build the 4x4 T matrix by putting the T_2x2 matrices on the diagonal
T_sectioned{1,1} = T_2x2;
T_sectioned{2,2} = T_2x2;
%rebuild the T matrix from the sectioned T matrix
T=cell2mat(T_sectioned);
%Find the 4 x 4 element stiffness matrix in global coordinates [KelGlobal]=[Ttranspose][KelLocal][T]
Kelglobal = transpose(T)*Kellocal*T;
%Identify the 2 dofs associated with node 1 of the element
dof_node1_index = find(dof_index(:,2)==node_1); %find all instances of node_1 id within the node column of the dof_index matrix
dof1_node1 = dof_index(dof_node1_index(1),1); %identify first dof associated with node 1
dof2_node1 = dof_index(dof_node1_index(2),1); %identify second dof associated with node 1
%Identify the 2 dofs associated with node 2 of the element
dof_node2_index = find(dof_index(:,2)==node_2); %find all instances of node_2 id within the node column of the dof_index matrix
dof1_node2 = dof_index(dof_node2_index(1),1); %identify first dof associated with node 2
dof2_node2 = dof_index(dof_node2_index(2),1); %identify second dof associated with node 2
%Partition [KelGlobal] into 4 2 x 2 matrices
Kelglobal_sectioned = mat2cell(Kelglobal,[2 2],[2 2]);
%Add element global stiffness matrix to global stiffness matrix in the
%correct locations
Kglobal_sectioned{node_1,node_1} = Kglobal_sectioned{node_1,node_1} + Kelglobal_sectioned{1,1};
Kglobal_sectioned{node_1,node_2} = Kglobal_sectioned{node_1,node_2} + Kelglobal_sectioned{1,2};
Kglobal_sectioned{node_2,node_1} = Kglobal_sectioned{node_2,node_1} + Kelglobal_sectioned{2,1};
Kglobal_sectioned{node_2,node_2} = Kglobal_sectioned{node_2,node_2} + Kelglobal_sectioned{2,2};

%initialize an element mass matrix
Melglobal = zeros(4,4);
%adding the element mass to the diagonal of the element mass matrix,
%corresponding to the 4 DOFs for n1 and n2 that the element connects
Melglobal(1,1) = m/2;
Melglobal(2,2) = m/2;
Melglobal(3,3) = m/2;
Melglobal(4,4) = m/2;
%Partition [Melglobal] into 4 2 x 2 matrices
Melglobal_sectioned = mat2cell(Melglobal,[2 2],[2 2]);
%adding the element mass matrix into the global mass matrix in correct
%location
Mglobal_sectioned{node_1,node_1} = Mglobal_sectioned{node_1,node_1} + Melglobal_sectioned{1,1};
Mglobal_sectioned{node_1,node_2} = Mglobal_sectioned{node_1,node_2} + Melglobal_sectioned{1,2};
Mglobal_sectioned{node_2,node_1} = Mglobal_sectioned{node_2,node_1} + Melglobal_sectioned{2,1};
Mglobal_sectioned{node_2,node_2} = Mglobal_sectioned{node_2,node_2} + Melglobal_sectioned{2,2};

%create temporary c vector
c_vectemp = zeros(1,n_dof);
%populate the c vector according to the element id to be able to calculate appropriate stresses for this element 
%populate the first trig term in the c vector
c_vectemp(1,dof1_node1)=-cos(theta);
%populate the second trig term in the c vector
c_vectemp(1,dof2_node1)=-sin(theta);
%populate the third trig term in the c vector
c_vectemp(1,dof1_node2)=cos(theta);
%populate the fourth trig term in the c vector
c_vectemp(1,dof2_node2)=sin(theta);
%write temporary c vector into overall c matrix
c(iel,1:n_dof) = c_vectemp;
end;
clear iel; %clears out the id variable so that it doesnt mess with indexing later

%Rebuild Kglobal from sectioned matrix
Kglobal = cell2mat(Kglobal_sectioned);

%Rebuild Mglobal from sectioned matrix
Mglobal = cell2mat(Mglobal_sectioned);

%WE NOW HAVE GLOBAL STIFFNESS AND MASS MATRICES OF THE STRUCTURE!! :)

%% Addition of the effect of non-structural masses

%read over the next two lines, the first which is not used for
%anything other than formatting and last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%Read in Number of concentrated non structural masses as an integer
NconcM = fscanf(fid,'%i',1);

%read over the next three lines, the first two which are not used for
%anything other than formatting and third last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%read in the properties for each concentrated mass and save in the array
%masses
concmasses=zeros(NconcM,2); %initialize the concmasses array, which has two columns for mass concmass(p) and the node it is attached to nodeconcmass(p)
for i=1:NconcM
 temp_ar=fscanf(fid,'%i %f %i',3); %reads in each concentrated mass into a temporary array consisting of concmass and nodeconcmass
 concmasses(temp_ar(1),1:2)=temp_ar(2:3); %saves the temporary array to the id'th row of the concmasses matrix
end
clear i; %clears out the id variable so that it doesnt mess with indexing later
%%
%Step 4
%Add the concentrated masses to the global mass matrix
for p=1:NconcM
%find the mass value 
Mp = concmasses(p,1);
%find the node the mass is located at
nodep = concmasses(p,2);
%find the associated degrees of freedom for the node
dof_node1_index = find(dof_index(:,2)==nodep); %find all instances of nodep id within the node column of the dof_index matrix
i1 = dof_index(dof_node1_index(1),1); %identify first dof associated with node
i2 = dof_index(dof_node1_index(2),1); %identify second dof associated with node
%add the concentrated masses to the appropriate DOFs
Mglobal(i1,i1)=Mglobal(i1,i1)+Mp;
Mglobal(i2,i2)=Mglobal(i2,i2)+Mp;
end

%% Forces and Boundary Conditions

%read over the next two lines, the first which is not used for
%anything other than formatting and last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%Read in Number of force cases as an integer
Nforcecases = fscanf(fid,'%i',1);

%read over the next three lines, the first two which are not used for
%anything other than formatting and third last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%read in the properties for each force case and save in the array
%forcecases
forcecases=zeros(Nforcecases,n_dof); %initialize the forcecases array
for i=1:Nforcecases
 temp_ar=fscanf(fid,'%f',n_dof+1); %reads in each force case into a temporary array consisting of its id and the vector of forces
 forcecases_transpose(temp_ar(1),1:n_dof)=temp_ar(2:n_dof+1); %saves the temporary array to the id'th row of the forcecases matrix
end
forcecases = forcecases_transpose'; %takes the transpose of forcecases_transpose so that each column in the forcecases matrix is a force vector

clear i; %clears out the id variable so that it doesnt mess with indexing later

%read over the next two lines, the first which is not used for
%anything other than formatting and last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%Read in Number of boundary condition sets as an integer
NBCsets = fscanf(fid,'%i',1);

%read over the next three lines, the first two which are not used for
%anything other than formatting and third last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%read in the properties for each concentrated mass and save in the array
%masses
BCdofID=zeros(NBCsets,n_dof); %initialize the BOCdofID array, which identifies if a DOF is constrained
for i=1:NBCsets
 temp_ar=fscanf(fid,'%f',n_dof+1); %reads in each constraint set into a temporary array consisting of its id and the vector of constraints
 BCdofID(temp_ar(1),1:n_dof)=temp_ar(2:n_dof+1); %saves the temporary array to the id'th column of the constraint set matrix
end
BCdofID = BCdofID'; %transpose the BCdofID array so the boundary conditions are column vectors within the array
clear i; %clears out the id variable so that it doesnt mess with indexing later

%% creating the NBCsets of global stiffness and mass matrices, each corresponding to a different boundary conditions (BC) set
%%
%Step 5
Kglobal_bcs = cell(NBCsets,1); %initialize a cell array to store the global stiffness matrix for each set of boundary conditions
Kglobal_backup = Kglobal; %initialize a backup global stiffness matrix 
Mglobal_bcs = cell(NBCsets,1); %initialize a cell array to store the global mass matrix for each set of boundary conditions
Mglobal_backup = Mglobal; %initialize a backup global mass matrix 
for ibc=1:NBCsets
 Kglobal = Kglobal_backup; %Restore the backup of the stiffness matrix so that each boundary condition can be applied on a virgin stiffness matrix
 Mglobal = Mglobal_backup; %Restore the backup of the mass matrix so that each boundary condition can be applied on a virgin mass matrix 
 temp_boundary=BCdofID(1:n_dof,ibc); %store the vector of the current boundary condition set
 for i=1:n_dof
     diag_elem=Kglobal(i,i); %store the diagonal element i,i of the global stiffness matrix 
     if temp_boundary(i)==0 
         Kglobal(i,1:n_dof)=0; %zero out the ith row of the stiffness matrix
         Kglobal(1:n_dof,i)=0; %zero out the ith column of the stiffness matrix
         Kglobal(i,i)=diag_elem; %restore the i,i diagonal element to the stiffness matrix
         Mglobal(i,1:n_dof)=0; %zero out the ith row of the mass matrix
         Mglobal(1:n_dof,i)=0; %zero out the ith column of the mass matrix
         Mglobal(i,i)=0.0000000001; %enter a small value of mass on the diagonal, leading to uncoupled natural frequencies of the dummy 1dof oscillators corresponding to constrained dof that will be high and out of range of interest (simple trick to deal with zero motion dofs)
     end
 end
 Kglobal_bcs{ibc}=Kglobal; %Store the global stiffness matrix corresponding to the ibc'th boundary condition set into the stiffness matrix cell array
 Mglobal_bcs{ibc}=Mglobal; %Store the global mass matrix corresponding to the ibc'th boundary condition set into the mass matrix cell array
end
clear i; %clears out the id variable so that it doesnt mess with indexing later
clear ibc; %clears out the id variable so that it doesnt mess with indexing later

%% Static Load cases

%read over the next two lines, the first which is not used for
%anything other than formatting and last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%Read in Number of static load cases sets as an integer
Nstatic = fscanf(fid,'%i',1);

%read over the next three lines, the first two which are not used for
%anything other than formatting and third last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%read in the properties for each concentrated mass and save in the array
%masses
static_load=zeros(Nstatic,2); %initialize the static_load array, which holds pointers to the force vector iforce and boundary condition case ibc for each static load case
for i=1:Nstatic
 temp_ar=fscanf(fid,'%i %i %i',3); %reads in each pointer into a temporary array 
 static_load(temp_ar(1),1:2)=temp_ar(2:3); %assigns the iforce and ibc from the temporary array to the id'th row of the static_load matrix
end
clear i; %clears out the id variable so that it doesnt mess with indexing later

%%
%Step 6

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

%% Solve the static problem [K]{u}={F}

%%
%Step 7
%find displacements
u = zeros(n_dof,Nstatic); %initialize displacement matrix to hold displacements for each static load case
stress_matrix = zeros(Nelements,Nstatic); %initialize stress matrix to hold stresses for each static load case and each element
for ist=1:Nstatic
    ibc = static_load(ist,2); %find the boundary condition set index ibc
    iforce = static_load(ist,1); %find the force set index iforce
    force = force_inputs(1:n_dof,iforce); %store the modified force vector associated to the index
    Kglobal_set = Kglobal_bcs{ibc}; %store the global stiffness matrix corresponding to BC set number ibc 
    u_case = lsqr(Kglobal_set,force); %Solve [K]{u} = {F} for {u}, lsqr function is used since K is sparse
    u(1:n_dof,Nstatic)=u_case; %store the displacement for this load case in the Nstatic column of the displacement matrix
    %find stresses
    stress_vec = zeros(1,Nelements); %initialize stress vector to capture stresses for each element
    for iel=1:Nelements
    c_element = c(iel,1:n_dof); %get the vector {c} for the element
    %Find the id_material for the element
    id_material = elements(iel,3);
    %Find material properties based on the id_material
    E=materials(id_material,1);
    %find ids of node 1 and node 2 that the element connects from elements
    %matrix
    node_1 = elements(iel,1);
    node_2 = elements(iel,2);
    %find x and y coordinates of node 1 and node 2 from the nodes matrix
    node_1x = nodes(node_1,1);
    node_1y = nodes(node_1,2);
    node_2x = nodes(node_2,1);
    node_2y = nodes(node_2,2);
    %calculate length of element 
    l = sqrt((node_2y-node_1y)^2+(node_2x-node_1x)^2);
    stress = (E/l)*c_element*u_case; %compute stress {c}*{u}
    stress_vec(iel)= stress; %save the stress for the iel element into the iel'th element of the stress_vec
    end
    stress_matrix(1:Nelements,ist) = stress_vec;
end

%% Static solution (Stresses in each element)
disp('The stresses in each element in PSI are displayed below in order of element number');
disp(stress_matrix);

%% Dynamic Cases

%read over the next two lines, the first which is not used for
%anything other than formatting and last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%Read in Number of elements as an integer
Ndynamic = fscanf(fid,'%i',1);

%read over the next three lines, the first two which are not used for
%anything other than formatting and third last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);
temp_str=fgetl(fid);

dynamic_cases=zeros(Ndynamic,2); %initialize the dynamic cases array, which has the 2 properties which we are tracking
for idyn=1:Ndynamic
 temp_ar=fscanf(fid,'%i %f %f',3); %reads in each element into a temporary array consisting of the boundary condition set to be used and the number of lowest natural frequencies we want to consider
 dynamic_cases(temp_ar(1),1:2)=temp_ar(2:3); %assigns the properties from the temporary array to the id'th row of the elements matrix
 ibc = dynamic_cases(idyn,1); %find the boundary condition set index ibc
 nlow = dynamic_cases(idyn,2); %find the number of lowest natural frequencies we want to concider, nlow
 Kglobal_set = Kglobal_bcs{ibc}; %store the global stiffness matrix corresponding to BC set number ibc 
 Mglobal_set = Mglobal_bcs{ibc}; %store the global mass matrix corresponding to BC set number ibc 
 [phi,lambda] = eig(Kglobal_set,Mglobal_set); %Solve the natural vibrations generalized eigenvalue problem lambda[M]{phi}=[K]{phi}
 lambdas_sorted = sort(lambda(lambda>0)); %Sort the lambda array from lowest to highest and eliminate zeros
 natural_freqs_sorted = (lambdas_sorted.^(1/2)).*(1/2*pi); %Convert the lambdas to natural frequencies in hertz
 natural_freqs = natural_freqs_sorted(1:nlow); %Chooes the nlow lowest natural frequencies and store in a vector
end
clear idyn; %clears out the id variable so that it doesnt mess with indexing later



