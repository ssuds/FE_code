function [design_variables, Nvars] = optinputs(file_name, Nelements)
    %Read input file for optimizer
    fid=fopen(file_name,'r'); 
    temp_str=fgetl(fid); %read over the first line, which is not used for anything other than formatting

    %read in number of design variable sets as an integer
    Nvars=fscanf(fid,'%i',1); 
    
    if Nvars == 0
        design_variables = []; %if there are no design variables, then return an empty design variables array which will tell the finite element code to run the nominal case 
    else
        %read over the next three lines, the first two which are not used for
        %anything other than formatting and third last being the one we want
        temp_str=fgetl(fid);
        temp_str=fgetl(fid);
        temp_str=fgetl(fid);

        %read in the design variables for each element and save in the array design_variables
        design_variables=zeros(Nelements,Nvars); %initialize the design variables array, which identifies if a particular characteristic of an element is a design variable
        for iel=1:Nelements
            temp_ar=fscanf(fid,'%f',Nvars+1); %reads in each elements design variables into a temporary array consisting of its id and the design variables
            design_variables(temp_ar(1),1:Nvars)=temp_ar(2:Nvars+1); %saves the temporary array to the id'th column of the design variable matrix
        end
    end
end