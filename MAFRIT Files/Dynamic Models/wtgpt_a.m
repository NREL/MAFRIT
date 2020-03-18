classdef wtgpt_a < Components
    % WTGPT_A represents the pitch controller module for Type 3 Wind Turbines
    % Based on the model description given here: https://www.wecc.org/Reliability/WECC-Second-Generation-Wind-Turbine-Models-012314.pdf    

    
    properties
        mvab;
        pref0;
        kiw;
        kpw;
        kic;
        kpc;
        kcc;
        tpi;
        pimax;
        pimin;
        piratmax;
        piratmin;
        wtgptStates
        wtgptOutputLocation;
        stWrite;
    end
    
    methods
        function wtgpt=wtgpt_a()
        end
        
        function wtgpt=AssignParams(wtgpt)
            wtgpt.mvab=wtgpt.compParams(1);
            wtgpt.kiw=wtgpt.compParams(2);
            wtgpt.kpw=wtgpt.compParams(3);
            wtgpt.kic=wtgpt.compParams(4);
            wtgpt.kpc=wtgpt.compParams(5);
            wtgpt.kcc=wtgpt.compParams(6);
            wtgpt.tpi=wtgpt.compParams(7);
            wtgpt.pimax=wtgpt.compParams(8);
            wtgpt.pimin=wtgpt.compParams(9);
            wtgpt.piratmax=wtgpt.compParams(10);
            wtgpt.piratmin=wtgpt.compParams(11);
            wtgpt.stWrite=fopen(wtgpt.wtgptOutputLocation,'w');
            if wtgpt.stWrite<0
                sprintf('error opening file')
            else
                fprintf(wtgpt.stWrite,'%8s %8s %8s \r\n',['State 1' 'State 2' 'State 3']);
            end  
        end
        %% Writing exciter states    
        function wtgpt=WriteData(wtgpt)
                fprintf(wtgpt.stWrite,'%8.4f %8.4f %8.4f \r\n',[wtgpt.wtgptStates(1) wtgpt.wtgptStates(2) wtgpt.wtgptStates(3)]);  
        end
        %%
        function [wtgpt,wt]=Init(wtgpt,pord,wref,theta0)
            wtgpt.pref0=pord;
            wt=wref;
            y7=theta0;
            y3=y7/2;
            y4=y3;     % Since y3+y4=y7 and independent assignment of y3 and y4 to obtain y7 is not possible equal values equal to half of y7 are assigned to y3 and y4.
            wtgpt.wtgptStates=[y3;y4;y7]; % 3 X 1 vector of initial conditions
        end
        %%
        function [wtgpt,theta]=newStates(wtgpt,tspan,pord,wref,wt)
            options=odeset('RelTol',1e-9);
            [T,X]=ode23t(@(t,x) wtgpt_dyn(t,x,pord,wref,wt,wtgpt.pref0,wtgpt.kiw,wtgpt.kpw,wtgpt.kic,wtgpt.kpc,wtgpt.kcc,wtgpt.tpi,wtgpt.pimax,wtgpt.pimin,wtgpt.piratmax,wtgpt.piratmin),tspan,wtgpt.wtgptStates,options);
            
            function dy=wtgpt_dyn(t1,y,pord,wref,wt,pref0,kiw,kpw,kic,kpc,kcc,tpi,pimax,pimin,piratmax,piratmin)
                
                dy=zeros(3,1);    
                y2=pord-pref0;
                y1=-wref+wt+kcc*y2;
                
                dy(1)=kiw*y1;
                if (y(1)>=pimax && dy(1)>0) || (y(1)<=pimin && dy(1)<0)
                    dy(1)=0;
                end
                y5=y(1)+kpw*y1;
                pi10=y5;
                if pi10>=pimax
                    pi10=pimax;
                elseif pi10<=pimin
                    pi10=pimin;
                end
                
                dy(2)=kic*y2;
                if (y(2)>=pimax && dy(2)>0) || (y(2)<=pimin && dy(2)<0)
                    dy(2)=0;
                end
                y6=y(2)+kpc*y2;
                pi20=y6;
                if pi20>=pimax
                    pi20=pimax;
                elseif pi20<=pimin
                    pi20=pimin;
                end
                
                y7=pi10+pi20;
                dy(3)=(y7-y(3))/tpi;
                if y(3)>=pimax && dy(3)>0
                    dy(3)=0;
                elseif y(3)<=pimin && dy(3)<0
                    dy(3)=0;
                elseif (y(3)<pimax && y(3)>pimin) && dy(3)>=piratmax
                    dy(3)=piratmax;
                elseif (y(3)<pimax && y(3)>pimin) && dy(3)<=piratmin
                    dy(3)=piratmin;
                end    
            end
            
            wtgpt.wtgptStates=X(end,:);   
            theta=wtgpt.wtgptStates(3);
        end
        
    end
    
end

