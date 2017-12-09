function [Nnodes,n_dof,nodes,Nmaterials,materials,Nelements,elements,NconcM,concmasses,Nforcecases,forcecases,NBCsets,BCdofID,Nstatic,static_load,Nmasscases,masscases,Ndynamic,dynamic_cases,sensitivityflag] = inputs(file_name)

%% Parse through input deck and store values into variables

%% Inputs for global stiffness and structural mass matrices

%Read input file, this is basically my punchcard!
fid=fopen(file_name,'r'); 
temp_str=fgetl(fid); %read over the first line, which is not used for anything other than formatting

%read in number of nodes as an integer
Nnodes=fscanf(fid,'%i',1); 

%Define degrees of freedom (2 times the number of nodes for our 2d anal)
n_dof = Nnodes*2;

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

%% Inputs for non-structural masses

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

%% Inputs for Forces and Boundary Conditions

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

%% Inputs for Static Load cases

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

%% Inputs for Mass Matrix Generation

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
for i=1:Nmasscases
    temp_str=fgetl(fid); %reads in each mass case line into a temporary string 
    temp_ar=sscanf(temp_str,'%f',NconcM+1); %breaks the temporary string into an array consisting of its id and all the associated concentrated mass ids
    masscases{i} = temp_ar(2:end); %places the mass case into the masscases cell array
end    

%% Inputs for Dynamic Cases

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

for idyn=1:Ndynamic
     temp_ar=fscanf(fid,'%i %f %f %f',4); %reads in each element into a temporary array consisting of the boundary condition set to be used, the mass case to be used and and the number of lowest natural frequencies we want to consider
     dynamic_cases(temp_ar(1),1:3)=temp_ar(2:4); %assigns the properties from the temporary array to the id'th row of the elements matrix
end

%% Read in sensitivity flag (Determines whether analytic sensitivites will be computed)

%read over the next two lines, the first which is not used for
%anything other than formatting and last being the one we want
temp_str=fgetl(fid);
temp_str=fgetl(fid);

%Read in sensitivity flag as an integer
sensitivityflag = fscanf(fid,'%i',1);

%Close the input deck
fclose(fid);

end
