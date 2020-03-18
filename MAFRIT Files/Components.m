classdef Components
    %COMPONENTS Parent Class for Generators and their Controllers
    %   COMPONENTS Class defines four properties that are common to any generator or its controller. These are:
    %       compParamLocation:      This variable stores the address of the text file that contains component parameters.
    %       compOutputLocation:     This variable stores the address of the text file where states of the components will be recorded at each iteration.
    %       compParams:             This variable stores a component's parameters.
    %       bus:                    This variable stores the bus # or ID at which a generator is connected. It also acts as the "primary key" for associating a generator with its controllers
    %   The COMPONENTS Class also includes two functions to read the
    %   parameters of components and write their states at the end of every
    %   simulation itertaion. These functions are:
    %       [comp,stRead]=ReadData(comp):   This function reads a component's parameters if the parameter file is successfully opened; otherwise it prints an error message.
    %       [comp,stWrite]=WriteData(comp): This function writes a component's parameters if the output file is successfully opened; otherwise it prints an error message.
    properties
         compParamLocation
         compOutputLocation
         compParams
         bus 
    end
    
    methods
        function comp=Components()
        end
        
        function [comp,stRead]=ReadData(comp)
            stRead=fopen(comp.compParamLocation,'r');
            if stRead<0
                sprintf('error opening file')
            else
                temp=cell2mat(textscan(stRead,'%f','delimiter',',','commentStyle','%'));
                comp.compParams=temp;
                fclose(stRead);
            end
        end     
    end
end


