classdef agcIsol < Components
    % AGC simple AGC model for an isolated system
  
    properties
        fs;
        tfilt;
        B;
        kp;
        ki;
        agcStates;
        agcOutputLocation;
        stWrite;
    end
    
    methods
        %% Constructor
        function agc=agcIsol()
        end
        %% Assigning AGC parameters to variables
        function agc=AssignParams(agc)
            agc.fs=agc.compParams(1);
            agc.tfilt=agc.compParams(2);
            agc.B=agc.compParams(3);
            agc.kp=agc.compParams(4);
            agc.ki=agc.compParams(5);
            agc.stWrite=fopen(agc.agcOutputLocation,'w');
            if agc.stWrite<0
                sprintf('error opening file')
            else
                fprintf(agc.stWrite,'%8s %8s \r\n',['State 1' 'State 2']);
            end           
        end
       %% Writing AGC states    
        function agc=WriteData(agc)
            fprintf(agc.stWrite,'%8.4f %8.4f \r\n',[agc.agcStates(1) agc.agcStates(2)]);     
        end
       %% Initilizing AGC states
        function [agc]=Init(agc)
            ace=0;
            agc.agcStates=[0;0]; % 2 X 1 vector of initial conditions
        end
        %% function is called at every simulation iteration to calculate new AGC states.
        function [agc,ace] = newStates(agc,tspan,fm)
            options=odeset('RelTol',1e-9);
            [T,X]=ode23t(@(t,x) agc_dyn(t,x,fm,agc.fs,agc.tfilt,agc.B,agc.ki),tspan,agc.agcStates,options);
            function dx = agc_dyn(t,x,fm,fs,tfilt,B,ki)
                dx=zeros(2,1); 
                ferr=-fm/60+fs;
                if tfilt==0
                    dx(1)=0;
                    perr=10*B*ferr;
                else
                    dx(1)=(ferr-x(1))/tfilt;
                    perr=10*B*x(1);
                end
                dx(2)=ki*perr;
            end
            agc.agcStates=X(end,:);
            perror=10*agc.B*agc.agcStates(1);
            ace=agc.kp*perror+agc.agcStates(2);
       end
    end   
end

