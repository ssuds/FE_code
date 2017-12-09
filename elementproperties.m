function [l,theta,E,A,V,rho,m,sigma_yt,sigma_yc] = elementproperties(n_dof,Nelements,elements,nodes,dof_index,materials)
% ELEMENTPROPERTIES  Computes commonly used properties for all elements in a truss structure.
%   

l = zeros(Nelements,1); %initialize empty length matrix to store lengths of each element
theta = zeros(Nelements,1); %initialize empty length matrix to store angle theta of each element
E = zeros(Nelements,1); %initialize empty matrix to store E of each element
sigma_yt = zeros(Nelements,1); %initialize empty matrix to store tensile yield stress allowable of each element
sigma_yc = zeros(Nelements,1); %initialize empty matrix to store compressive yield stress allowable of each element
A = zeros(Nelements,1); %intiialize empty matrix to store A of each element
V = zeros(Nelements,1); %intiialize empty matrix to store V of each element
rho = zeros(Nelements,1); %initialize empty length matrix to store lengths of each element
m = zeros(Nelements,1); %intiialize empty matrix to store mass of each element
for iel=1:Nelements
    %find ids of node 1 and node 2 that the element connects from elements matrix
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
    %calculate the angle the line from node 1 to node 2 creates with the global x axis
    theta_el=atan2(node_2y-node_1y,node_2x-node_1x);
    theta(iel,1) = theta_el; %store the angle of this element in the theta vector
    %Find the id_material for the element
    id_material = elements(iel,3);
    %Find material properties based on the id_material
    E_el=materials(id_material,1);
    E(iel,1) = E_el; %store E of this element in the E vector
    sigma_yt_el=materials(id_material,2);
    sigma_yt(iel,1) = sigma_yt_el; %store tensile yield stress allowable of this element in the sigma_yt vector
    sigma_yc_el=materials(id_material,3);
    sigma_yc(iel,1) = sigma_yc_el; %store compressive yield stress allowable of this element in the sigma_yc vector
    rho_el=materials(id_material,4);
    rho(iel,1) = rho_el; %store the density rho of this element 
    %Find the cross sectional area of the element A-iel
    A_el = elements(iel,4);
    A(iel,1) = A_el;
    %Find the volume of the element A-iel * L-iel (note: the –iel is an index
    %not a subtraction)
    V_el = A_el*l_el;
    V(iel,1) = V_el;
    %Find the mass of the element m-iel=elementMaterialDensity*elementVolume
    m_el=rho(iel,1)*V(iel,1);
    m(iel,1) = m_el;
end


end