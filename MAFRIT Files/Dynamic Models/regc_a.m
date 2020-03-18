classdef regc_a < Components
    % regc_A represents the electrical control module for Type 3 and 4 Wind Turbines and Central PV Stations.
    % Based on the model description given here: https://www.wecc.org/Reliability/WECC-Second-Generation-Wind-Turbine-Models-012314.pdf    
 
    
    properties
        mvab;
        tfltr;
        lvpl1;
        zerox;
        brkpt;
        lvplsw;
        rrpwr;
        tg;
        volim;
        iolim;
        khv;
        ivpnt0;
        ivpnt1;
        iqrmax;
        iqrmin;
        xd11;
        regcStates;
        regcOutputLocation;
        stWrite;
    end
    
    methods
        function regc=regc_a()
        end
        
        function regc=AssignParams(regc)
            regc.mvab = regc.compParams(1);
            regc.tfltr = regc.compParams(2);
            regc.lvpl1 = regc.compParams(3);
            regc.zerox = regc.compParams(4);
            regc.brkpt = regc.compParams(5);
            regc.lvplsw = regc.compParams(6);
            regc.rrpwr = regc.compParams(7);
            regc.tg = regc.compParams(8);
            regc.volim = regc.compParams(9);
            regc.iolim = regc.compParams(10);
            regc.khv = regc.compParams(11);
            regc.ivpnt0 = regc.compParams(12);
            regc.ivpnt1 = regc.compParams(13);
            regc.iqrmax = regc.compParams(14);
            regc.iqrmin = regc.compParams(15);
            regc.xd11 = regc.compParams(16);
            regc.stWrite=fopen(regc.regcOutputLocation,'w');
            if regc.stWrite<0
                sprintf('error opening file')
            else
                fprintf(regc.stWrite,'%12s %12s %12s \r\n',['State 1' 'State 2' 'State 3']);
            end  
        end
        %% Writing exciter states    
        function regc=WriteData(regc)
                fprintf(regc.stWrite,'%8.4f %8.4f %8.4f \r\n',[regc.regcStates(1) regc.regcStates(2) regc.regcStates(3)]);  
        end
        %%
        function [regc,ipcmd,iqcmd]=Init(regc,inet,vt,qgen)
            vterm=abs(vt);
            iq2=-imag(conj(inet)*vt);
            ip1=real(conj(inet)*vt);
            v=vterm;
            ipcmd=ip1/vterm;
            iqcmd=-iq2/vterm;
            if qgen>=0
                regc.iqrmin=-9999;  % deactivate lower reactive current rate limit when initial qgen >0
            elseif qgen<0
                regc.iqrmax=9999;   % deactivate upper reactive current rate limit when initial qgen <0
            end
            regc.regcStates=[-iqcmd;v;ipcmd]; % 3 X 1 vector of initial conditions
        end
        %%
        function [regc,inet]=newStates(regc,tspan,vt,ipcmd,iqcmd)
            options=odeset('RelTol',1e-9);
            [T,X]=ode23t(@(t,x) regc_dyn(t,x,ipcmd,iqcmd,vt,regc.tfltr,regc.lvpl1,regc.zerox,regc.brkpt,regc.lvplsw,regc.rrpwr,regc.tg,regc.iqrmax,regc.iqrmin),tspan,regc.regcStates,options);
            
            function dy=regc_dyn(t1,y,ipcmd,iqcmd,vt,tfltr,lvpl1,zerox,brkpt,lvplsw,rrpwr,tg,iqrmax,iqrmin)
                vterm=abs(vt);
                dy=zeros(3,1); 
                dy(1)=(-iqcmd-y(1))/tg;
                if dy(1)>=iqrmax
                    dy(1)=iqrmax;
                elseif dy(1)<=iqrmin
                    dy(1)=iqrmin;
                end  

                if lvplsw==0
                    lvpl=99;
                    dy(3)=0;
                elseif lvplsw==1
                    dy(2)=(vterm-y(2))/tfltr;
                    if y(2)<=zerox
                        lvpl=0;
                    elseif y(2)>zerox && y(2)<brkpt
                        lvpl=lvpl1*(y(2)-zerox)/(brkpt-zerox);
                    else
                        lvpl=99;
                    end
                end
                
                dy(3)=(ipcmd-y(3))/tg;
                if lvplsw==1 && y(3)>=lvpl && dy(3)>0
                    dy(3)=0;
                elseif lvplsw==1 && y(3)<lvpl && dy(3)>rrpwr
                    dy(3)=rrpwr;
                end
                    
            end
            regc.regcStates=X(end,:);          
            iq=regc.regcStates(1);
            v=regc.regcStates(2);
            ip=regc.regcStates(3);
            iq2=abs(vt)*iq;
            v2=abs(vt)-regc.volim;
            if abs(vt)<=regc.volim
                v3=0;
            else
                v3=regc.khv*v2;
            end
            if (iq2-v3)>regc.iolim
                iq1=iq2-v3;
            else
                iq1=regc.iolim;
            end
            
            if abs(vt)<=regc.ivpnt0
                gain=0;
            elseif abs(vt)>regc.ivpnt0 && abs(vt)<=regc.ivpnt1
                gain=(v-regc.ivpnt0)/(regc.ivpnt1-regc.ivpnt0);
            else
                gain=1;
            end
            ip1=ip*gain*abs(vt);
            inet=conj((ip1-1i*iq1)/vt); 
        end
        
    end
    
end

