classdef unitController < Components
    % unitController simple unit controller model for an isolated system to
    % allocate system-wide ACE to each unit
  
    properties
        allocFactor;
        rmax;
        rmin;
        kpref;
        uCntStates;
        uCntOutputLocation;
        stWrite;
    end
    
    methods
        %% Constructor
        function uCnt=unitController()
        end
        %% Assigning uCnt parameters to variables
        function uCnt=AssignParams(uCnt)
            uCnt.allocFactor=uCnt.compParams(1);
            uCnt.rmax=uCnt.compParams(2);
            uCnt.rmin=uCnt.compParams(3);
            uCnt.kpref=uCnt.compParams(4);
            uCnt.stWrite=fopen(uCnt.uCntOutputLocation,'w');
            if uCnt.stWrite<0
                sprintf('error opening file')
            else
                fprintf(uCnt.stWrite,'%8s \r\n',['State 1']);
            end           
        end
       %% Writing uCnt states    
        function uCnt=WriteData(uCnt)
            fprintf(uCnt.stWrite,'%8.4f \r\n',uCnt.uCntStates(1));     
        end
       %% Initilizing uCnt states
        function uCnt=Init(uCnt)
            uCnt.uCntStates=0; 
        end
        %% function is called at every simulation iteration to calculate new uCnt states.
        function [uCnt,dpref] = newStates(uCnt,tspan,ace,Sgov)
            options=odeset('RelTol',1e-9);
            [T,X]=ode23t(@(t,x) uCnt_dyn(t,x,ace,Sgov,uCnt.allocFactor,uCnt.rmax,uCnt.rmin),tspan,uCnt.uCntStates,options);
            function dx = uCnt_dyn(t,x,ace,Sgov,allocFactor,rmax,rmin)
                acei=allocFactor*ace;
                aceipu=acei/Sgov;
                if abs(aceipu)>=1
                   aceipu=aceipu/abs(aceipu);
                end
                dx(1)=aceipu-x(1);
                if dx(1)>rmax
                    dx(1)=rmax;
                elseif dx(1)<=rmin
                    dx(1)=rmin;
                end
            end
            uCnt.uCntStates=X(end,:);
            dpref=uCnt.uCntStates(1)*uCnt.kpref;
       end
    end   
end

