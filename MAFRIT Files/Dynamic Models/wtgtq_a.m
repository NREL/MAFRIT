classdef wtgtq_a < Components
    % WTGTQ_A represents the torque controller module for Type 3 Wind
    % Turbines
    % Based on the model description given here: https://www.wecc.org/Reliability/WECC-Second-Generation-Wind-Turbine-Models-012314.pdf    
    
    properties
        mvab;
        kip;
        kpp;
        tp;
        twref;
        temax;
        temin;
        p1;
        spd1;
        p2;
        spd2;
        p3;
        spd3;
        p4;
        spd4;
        tflag;
        wtgtqStates
        wtgtqOutputLocation;
        stWrite;
    end
    
    methods
        function wtgtq=wtgtq_a()
        end
        
        function wtgtq=AssignParams(wtgtq)
            wtgtq.mvab=wtgtq.compParams(1);
            wtgtq.kip=wtgtq.compParams(2);
            wtgtq.kpp=wtgtq.compParams(3);
            wtgtq.tp=wtgtq.compParams(4);
            wtgtq.twref=wtgtq.compParams(5);
            wtgtq.temax=wtgtq.compParams(6);
            wtgtq.temin=wtgtq.compParams(7);
            wtgtq.p1=wtgtq.compParams(8);
            wtgtq.spd1=wtgtq.compParams(9);
            wtgtq.p2=wtgtq.compParams(10);
            wtgtq.spd2=wtgtq.compParams(11);
            wtgtq.p3=wtgtq.compParams(12);
            wtgtq.spd3=wtgtq.compParams(13);
            wtgtq.p4=wtgtq.compParams(14);
            wtgtq.spd4=wtgtq.compParams(15);
            wtgtq.tflag=wtgtq.compParams(16);
            wtgtq.stWrite=fopen(wtgtq.wtgtqOutputLocation,'w');
            if wtgtq.stWrite<0
                sprintf('error opening file')
            else
                fprintf(wtgtq.stWrite,'%8s %8s %8s \r\n',['State 1' 'State 2' 'State 3']);
            end  
        end
        %% Writing exciter states    
        function wtgtq=WriteData(wtgtq)
                fprintf(wtgtq.stWrite,'%8.4f %8.4f %8.4f \r\n',[wtgtq.wtgtqStates(1) wtgtq.wtgtqStates(2) wtgtq.wtgtqStates(3)]);  
        end
        %%
        function [wtgtq,wg,pref0]=Init(wtgtq,pref,pe)
            pe1=pe;
            pr2=pref;
            if wtgtq.tflag==1
                pref0=pe;
            else
                pref0=0;      % implies pref0 is not used during control
            end
            
            if pe>=0 && pe<wtgtq.p1
                wref=wtgtq.spd1;
            elseif pe>=wtgtq.p1 && pe<wtgtq.p2
                wref=wtgtq.spd1+(pe-wtgtq.p1)*(wtgtq.spd2-wtgtq.spd1)/(wtgtq.p2-wtgtq.p1);
            elseif pe>=wtgtq.p2 && pe<wtgtq.p3
                wref=wtgtq.spd2+(pe-wtgtq.p2)*(wtgtq.spd3-wtgtq.spd2)/(wtgtq.p3-wtgtq.p2);
            elseif pe>=wtgtq.p3 && pe<wtgtq.p4
                wref=wtgtq.spd3+(pe-wtgtq.p3)*(wtgtq.spd4-wtgtq.spd3)/(wtgtq.p4-wtgtq.p3);
            else
                wref=wtgtq.spd4;
            end
            wg=wref;      % wg is set equal to wref irrespective of the value of Tflag. This is done because wref is used as an input for wtgpt model and it expects an input whether Tflag is 0 or 1.
            
            wtgtq.wtgtqStates=[pe1;wref;pr2]; % 3 X 1 vector of initial conditions
        end
        %%
        function [wtgtq,pref]=newStates(wtgtq,tspan,pe,wg,voltage_dip,pref0)
            options=odeset('RelTol',1e-9);
            [T,X]=ode23t(@(t,x) wtgtq_dyn(t,x,pe,wg,voltage_dip,pref0,wtgtq.kip,wtgtq.tp,wtgtq.twref,wtgtq.temax,wtgtq.temin,wtgtq.p1,wtgtq.spd1,wtgtq.p2,wtgtq.spd2,wtgtq.p3,wtgtq.spd3,wtgtq.p4,wtgtq.spd4,wtgtq.tflag),tspan,wtgtq.wtgtqStates,options);
            
            function dy=wtgtq_dyn(t1,y,pe,wg,voltage_dip,pref0,kip,tp,twref,temax,temin,p1,spd1,p2,spd2,p3,spd3,p4,spd4,tflag)
                 dy=zeros(3,1);    
                
                dy(1)=(pe-y(1))/tp;
                p0=(pref0-y(1))/wg;
                
                if y(1)>=0 && y(1)<p1
                    w=spd1;
                elseif y(1)>=p1 && y(1)<p2
                    w=spd1+(y(1)-p1)*(spd2-spd1)/(p2-p1);
                elseif y(1)>=p2 && y(1)<p3
                    w=spd2+(y(1)-p2)*(spd3-spd2)/(p3-p2);
                elseif y(1)>=p3 && y(1)<p4
                    w=spd3+(y(1)-p3)*(spd4-spd3)/(p4-p3); 
                else
                    w=spd4;
                end
                
                dy(2)=(w-y(2))/twref;
                
                w0=wg-y(2);
                
                if tflag==1
                    pr1=p0;
                elseif tflag==0
                    pr1=w0;
                end
                
                dy(3)=kip*pr1;
                if (y(3)>=temax && dy(3)>0) || (y(3)<=temin && dy(3)<0)
                    dy(3)=0;
                end
                
                if voltage_dip==1
                    dy(3)=0;
                end              
            end
            wtgtq.wtgtqStates=X(end,:); 
                p0_1=(pref0-wtgtq.wtgtqStates(1))/wg;               
                w0_1=wg-wtgtq.wtgtqStates(2);
                
                if wtgtq.tflag==1
                    pr1_1=p0_1;
                elseif wtgtq.tflag==0
                    pr1_1=w0_1;
                end
            pr3=wtgtq.wtgtqStates(3)+pr1_1*wtgtq.kpp;
            pref=pr3;
            if pr3>=wtgtq.temax
                pref=wtgtq.temax;
            elseif pr3<=wtgtq.temin
                pref=wtgtq.temin;
            end
        end
        
    end
    
end

