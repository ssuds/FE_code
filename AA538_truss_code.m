%% Truss FE code for 2D structures - AA 538 Introduction to Structural Optimization - Shreyas Sudhakar

%Clearing all previous data in workspace
clear; 

%% Contents of input deck
type input.txt;

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
% Step 1: Define the degrees of freedom and associate each dof with a node and a direction

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
% Step 2: Allocate memory for the global stiffness and mass matrices for the structure: [Kglobal],[Mglobal], initiate these matrices by setting all their elements to zero.

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

l = zeros(Nelements,1); %initialize empty length matrix to store lengths of each element
E = zeros(Nelements,1); %initialize empty matrix to store E of each element
A = zeros(Nelements,1); %intiialize empty matrix to store A of each element
rho = zeros(Nelements,1); %initialize empty length matrix to store lengths of each element
T_elements = cell(1,Nelements); %initialize cell matrix to store transpose matrix for each element
Kellocal_elements = cell(1,Nelements); %initialize cell matrix to store local stiffness matrix matrix for each element
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
l_el = sqrt((node_2y-node_1y)^2+(node_2x-node_1x)^2);
l(iel,1) = l_el; %store length of this element in the l vector
%calculate the angle the line from node 1 to node 2 creates with the
%global x axis
theta=atan2(node_2y-node_1y,node_2x-node_1x);
%Find the id_material for the element
id_material = elements(iel,3);
%Find material properties based on the id_material
E_el=materials(id_material,1);
E(iel,1) = E_el; %store E of this element in the E vector
rho_el=materials(id_material,4);
rho(iel,1) = rho_el; %store the density rho of this element 
%Find the cross sectional area of the element A-iel
A_el = elements(iel,4);
A(iel,1) = A_el;
%Find the volume of the element A-iel * L-iel (note: the –iel is an index
%not a subtraction)
V = A_el*l_el;
%Find the mass of the element m-iel=elementMaterialDensity*elementVolume
m=rho_el*V;
%generate the local stiffness matrix for the element
%note that A is multiplied in after the fact to streamline the code when
%optimization is in play

Kellocal = zeros(4,4);
Kellocal(1,1) = 1;
Kellocal(3,1) = -1;
Kellocal(1,3) = -1;
Kellocal(3,3) = 1;
Kellocal = (E_el/l_el).*Kellocal;
Kellocal = A_el.*Kellocal;
Kellocal_elements{iel} = Kellocal; %store the local stiffness matrix in a cell array
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
%store T matrix for each element
T_elements{iel} = T;
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
Melglobal = diag([1 1 1 1]);
Melglobal = m/2.*Melglobal;
%Partition [Melglobal] into 4 2 x 2 matrices
Melglobal_sectioned = mat2cell(Melglobal,[2 2],[2 2]);
%adding the element mass matrix into the global mass matrix in correct location
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

%% Forces and Boundary Conditions

%read over the next lines, which is not used for
%anything other than formatting
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
 forcecases_transpose(i,1:n_dof)=temp_ar(2:n_dof+1); %saves the temporary array to the id'th row of the forcecases matrix
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
%Step 6: create the force input vectors for each static load case

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
%Step 7: Solve for displacements and stresses in all static load cases
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
    E_el=materials(id_material,1);
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
    l_el = sqrt((node_2y-node_1y)^2+(node_2x-node_1x)^2);
    stress = (E_el/l_el)*c_element*u_case; %compute stress {c}*{u}
    stress_vec(iel) = stress; %save the stress for the iel element into the iel'th element of the stress_vec
    end
    stress_matrix(1:Nelements,ist) = stress_vec;
end

%% Static solution (Stresses in each element)
disp('The stresses in each element in PSI are displayed below in order of element number (vertically) and load case (horizontally)');
disp(stress_matrix);

%% HW2 part a) Add mass matrix generation

%read over the next two lines, the first which is not used for
%anything other than formatting and last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%Read in Number of mass case sets as an integer
Nmasscases = fscanf(fid,'%i',1);

%read over the next three lines, the first two which are not used for
%anything other than formatting and third last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%read in the concentrated masses in each mass case and save in the array
masscases={1,1:Nmasscases}; %initialize the masscases cell array
Mglobal_masscases = cell(Nmasscases,1); %initialize a cell array to store the global mass matrix for each set of concentrated mass cases
Mglobal_backup = Mglobal; %initialize a backup global mass matrix 
for i=1:Nmasscases
 temp_str=fgetl(fid); %reads in each mass case line into a temporary string 
 temp_ar=sscanf(temp_str,'%f',NconcM+1); %breaks the temporary string into an array consisting of its id and all the associated concentrated mass ids
 masscases{i} = temp_ar(2:end); %places the mass case into the masscases cell array
    for j=1:length(masscases{i})
        Mglobal = Mglobal_backup; %Restore the backup of the mass matrix so that each boundary condition can be applied on a virgin mass matrix 
        %find the mass value 
        Mj = concmasses(masscases{i}(j),1);
        %find the node the mass is located at
        nodej = concmasses(masscases{i}(j),2);
        %find the associated degrees of freedom for the node
        dof_node1_index = find(dof_index(:,2)==nodej); %find all instances of nodej id within the node column of the dof_index matrix
        i1 = dof_index(dof_node1_index(1),1); %identify first dof associated with node
        i2 = dof_index(dof_node1_index(2),1); %identify second dof associated with node
        %add the concentrated masses to the appropriate DOFs
        Mglobal(i1,i1)=Mglobal(i1,i1)+Mj;
        Mglobal(i2,i2)=Mglobal(i2,i2)+Mj;
    end 
 Mglobal_masscases{i}=Mglobal; %Store the global mass matrix corresponding to the i'th mass case set into the mass case matrix cell array
 end
clear i; %clears out the id variable so that it doesnt mess with indexing later
clear j; %clears out the id variable
%% HW2 part b) Define the number of natural frequency / mode shape cases
%% Dynamic Cases

%read over the next two lines, the first which is not used for
%anything other than formatting and last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%Read in Number of dynamic cases as an integer
Ndynamic = fscanf(fid,'%i',1);

%read over the next three lines, the first two which are not used for
%anything other than formatting and third last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);
temp_str=fgetl(fid);

dynamic_cases=zeros(Ndynamic,2); %initialize the dynamic cases array, which has the 2 properties which we are tracking
natural_freqs = cell(1,Ndynamic); %initialize a cell array to store the natural frequencies found for each dynamic case
phi = cell(1,Ndynamic); %initialize a cell array to store the eigenvectors phi found for each dynamic case
lambda = cell(1,Ndynamic); %initialize a cell array to store the eigenvalues lambda found for each dynamic case

for idyn=1:Ndynamic
 temp_ar=fscanf(fid,'%i %f %f %f',4); %reads in each element into a temporary array consisting of the boundary condition set to be used, the mass case to be used and and the number of lowest natural frequencies we want to consider
 dynamic_cases(temp_ar(1),1:3)=temp_ar(2:4); %assigns the properties from the temporary array to the id'th row of the elements matrix
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

%% HW2 part c) Analytic sensitivities with respect to cross sectional areas - static cases

%read over the next two lines, the first which is not used for
%anything other than formatting and last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%Read in sensitivity flag as an integer
sensitivityflag = fscanf(fid,'%i',1);

if sensitivityflag == 1 %Analytic sensitivities for will only be computed if the sensitivity flag is 1 within the input deck
    
    %initialize dK_dA cell matrix, which will hold the d[K]/dA used in computing analytic sensitivities for each element
    dK_dA = cell(Nelements,1);

    clear temp_ar;
    temp_ar(1:(n_dof/2))=2; 
    for iel = 1:Nelements
        dK_dAglobal = zeros(n_dof,n_dof); %initialize the global d[K]/dA matrix which is used for analytic sensitivities
        dK_dAglobal_sectioned = mat2cell(dK_dAglobal,temp_ar,temp_ar); %Section dK/dA matrix into 2x2 sections
        %build the element contribution to the global dK/dA matrix
        dK_dAlocal = zeros(4,4);
        dK_dAlocal(1,1) = 1;
        dK_dAlocal(3,1) = -1;
        dK_dAlocal(1,3) = -1;
        dK_dAlocal(3,3) = 1;
        dK_dAlocal = (E(iel,1)/l(iel,1)).*Kellocal_elements{iel};
        dK_dAelglobal = transpose(T_elements{iel})*dK_dAlocal*T_elements{iel}; %transform the element dK/dA matrix to the global coordinate system
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
        for iel=1:Nelements
            pseudo_load_loadcase(1:n_dof,iel) = dF_dA-(dK_dA{iel}*u_case); %compute pseudo load for each element
            du_dA_loadcase(1:n_dof,iel) = lsqr(Kglobal_set,pseudo_load_loadcase(1:n_dof,iel)); %Solve [K]{u} = {F} for {u}, lsqr function is used since K is sparse
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
        dM_dAelglobal = ((rho(iel,1)*l(iel,1))/2).*Melglobal; %Generate the element contribution to the global dM/dA matrix
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
 end
%% HW2 Part e) Truss member buckling

g_buckling = zeros(Nelements,Nstatic); %initialize matrix to hold buckling constraint values for all elements and static load cases
for ist = 1:Nstatic
    g_buckling_case = zeros(Nelements,1); %initialize vector to hold buckling constraint values for each element
    for iel = 1:Nelements
      g_buckling_el = -1-stress_matrix(iel,ist)*((4*l(iel,1)^2)/(pi*E(iel,1)*A(iel,1))); %compute buckling criterion (if this is less than or equal to zero, buckling occurs)
      g_buckling_case(iel,1) = g_buckling_el; %store the buckling constraint value for this element into the larger vector
    end
    
    dgbuckling_dAj = zeros(Nelements); %initialize array to store sensitivities of buckling constraints in all elements with respect to areas of all elements
   
if sensitivityflag == 1 %Analytic sensitivities will only be computed if the sensitivity flag is 1 within the input deck 
    %compute sensitivity of buckling constraints 
    for ieli = 1:Nelements
        for ielj = 1:Nelements
            if ielj == ieli %the second term in the sensitivity equation appears only when calculating the sensitivity of a buckling constraint in an element wrt the area of the same element
                dgbucklingi_dAj_ielj = (dsigma_dA{ist}(1,ielj)*((4*l(ieli,1)^2)/(pi*E(ieli,1)*A(ieli,1))))+(stress_matrix(ieli,ist)*((4*l(ieli,1)^2)/(pi*E(ieli,1)*(A(ielj,1)^2))));
            else
                dgbucklingi_dAj_ielj = dsigma_dA{ist}(1,ielj)*((4*l(ieli,1)^2)/(pi*E(ieli,1)*A(ieli,1)));
            end
            dgbuckling_dAj_ieli(ielj) = dgbucklingi_dAj_ielj; %store buckling constraints in one element with respect to areas of all elements
        end 
        dgbuckling_dAj(1:Nelements,ieli) = dgbuckling_dAj_ieli; %store buckling constraints in all elements with respect to areas of all elements
    end
end
            
    g_buckling(1:Nelements,ist) = g_buckling_case; %store the buckling constraint values for this static load case into the larger matrix
end
