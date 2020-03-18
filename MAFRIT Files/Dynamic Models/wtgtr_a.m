classdef wtgtr_a < Components
    % WTGTR_A represents the drive train model for Type 3 Wind Turbines.
    % The WTGAR_A Model which represents the wind turbine aerodynamics is
    % also combined with the wtgt model as there are no dynamic states in
    % this model and its output is fed only to the wtgt model.
    % Based on the model description given here: https://www.wecc.org/Reliability/WECC-Second-Generation-Wind-Turbine-Models-012314.pdf    

    
    properties
        mvab;
        ht;
        hg;
        dshaft;
        kshaft;
        ka;
        theta0;
        pm0;
        w0;
        wtgtrStates;
        wtgtrOutputLocation;
        stWrite;
    end
    
    methods
        function wtgtr=wtgtr_a()
        end
        
        function wtgtr=AssignParams(wtgtr)
            wtgtr.mvab=wtgtr.compParams(1);
            wtgtr.ht=wtgtr.compParams(2);
            wtgtr.hg=wtgtr.compParams(3);
            wtgtr.dshaft=wtgtr.compParams(4);
            wtgtr.kshaft=wtgtr.compParams(5);
            wtgtr.ka=wtgtr.compParams(6);
            wtgtr.theta0=wtgtr.compParams(7);
            wtgtr.w0=wtgtr.compParams(8);
            wtgtr.stWrite=fopen(wtgtr.wtgtrOutputLocation,'w');
            if wtgtr.stWrite<0
                sprintf('error opening file')
            elseif wtgtr.kshaft>0
                fprintf(wtgtr.stWrite,'%8s %8s %8s \r\n',['State 1' 'State 2' 'State 3']);
            elseif wtgtr.kshaft<=0
                fprintf(wtgtr.stWrite,'%8s \r\n',['State 1']);
            end  
        end
        %% Writing exciter states    
        function wtgtr=WriteData(wtgtr)
                if wtgtr.kshaft>0
                    fprintf(wtgtr.stWrite,'%8.6f %8.6f %8.6f \r\n',[wtgtr.wtgtrStates(1) wtgtr.wtgtrStates(2) wtgtr.wtgtrStates(3)]); 
                else
                    fprintf(wtgtr.stWrite,'%8.6f \r\n',[wtgtr.wtgtrStates(1)]); 
                end
        end
        %%
        function [wtgtr,theta0]=Init(wtgtr,pe,wref)
            if wtgtr.kshaft>0   % 2 Mass Model
                wtgtr.pm0=pe;    % Under the assumption that at the steady state theta=theta0
                wtgtr.w0=wref;        
                w1=0;
                w2=w1;
                T3=pe/wref;
                wtgtr.wtgtrStates=[w1;w2;T3/wtgtr.kshaft]; % 3 X 1 vector of initial conditions
                theta0=wtgtr.theta0;
            else
                wtgtr.w0=wref;
                wtgtr.pm0=pe;
                wtgtr.wtgtrStates=0; % 1 Mass Model
                theta0=wtgtr.theta0;
            end
        end
        %%
        function [wtgtr,wt,wg]=newStates(wtgtr,tspan,pe,theta)
            options=odeset('RelTol',1e-9);
            pm=wtgtr.pm0-wtgtr.ka*theta*(theta-wtgtr.theta0);
            [T,X]=ode23t(@(t,x) wtgtr_dyn(t,x,pe,wtgtr.ht,wtgtr.hg,wtgtr.dshaft,wtgtr.kshaft,wtgtr.w0),tspan,wtgtr.wtgtrStates,options);
            
            function dy=wtgtr_dyn(t1,y,pe,ht,hg,dshaft,kshaft,w0)
                if kshaft>0
                    dy=zeros(3,1); 
                    T1=(y(1)-y(2))*dshaft;
                    wt1=w0+y(1);
                    wg1=w0+y(2);
                    dy(1)=(pm/wt1-y(3)*kshaft-T1)/(2*ht);
                    if hg==0
                        dy(2)=0;
                    else
                        dy(2)=(-pe/wg1+T1+y(3)*kshaft)/(2*hg);
                    end
                    dy(3)=(y(1)-y(2));
                else
                    dy(1)=(1/(2*ht))*(pm-pe-y(1)*dshaft)/(w0+y(1));
                end
            end
            wtgtr.wtgtrStates=X(end,:);
            if wtgtr.kshaft>0
                wt=wtgtr.w0+wtgtr.wtgtrStates(1);
                wg=wtgtr.w0+wtgtr.wtgtrStates(2);
            else
                wt=wtgtr.w0+wtgtr.wtgtrStates(1);
                wg=wt;
            end
        end
        
    end
    
end

