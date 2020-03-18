classdef reec_b < Components
    % reec_b represents the electrical control module for Central PV Stations and battery storage plants.
    % Based on the model description given here: https://www.wecc.org/Reliability/WECC-Second-Generation-Wind-Turbine-Models-012314.pdf    
 
    properties
        mvab;
        vdip;
        vup;
        trv;
        dbd1;
        dbd2;
        kqv;
        iqhl;
        iqll;
        vref0;
        pfaref;  %this parameter is not input, but is assigned during initialization
        tp;
        qmax;
        qmin;
        vmax;
        vmin;
        kqp;
        kqi;
        kvp;
        kvi;
        tiq;
        dpmax;
        dpmin;
        pmax;
        pmin;
        imax;
        tpord;
        pfflag;
        vflag;
        qflag;
        pqflag;
        reecStates;
        reecOutputLocation;
        stWrite;
    end
    
    methods
        function reec=reec_b()
        end
        
        function reec=AssignParams(reec)
            reec.mvab = reec.compParams(1);
            reec.vdip = reec.compParams(2);
            reec.vup = reec.compParams(3);
            reec.trv = reec.compParams(4);
            reec.dbd1 = reec.compParams(5);
            reec.dbd2 = reec.compParams(6);
            reec.kqv = reec.compParams(7);
            reec.iqhl = reec.compParams(8);
            reec.iqll = reec.compParams(9);
            reec.vref0 = reec.compParams(10);
            reec.tp = reec.compParams(11);
            reec.qmax = reec.compParams(12);
            reec.qmin = reec.compParams(13);
            reec.vmax = reec.compParams(14);
            reec.vmin = reec.compParams(15);
            reec.kqp = reec.compParams(16);
            reec.kqi = reec.compParams(17);
            reec.kvp = reec.compParams(18);
            reec.kvi = reec.compParams(19);
            reec.tiq = reec.compParams(20);
            reec.dpmax = reec.compParams(21);
            reec.dpmin = reec.compParams(22);
            reec.pmax = reec.compParams(23);
            reec.pmin = reec.compParams(24);
            reec.imax = reec.compParams(25);
            reec.tpord = reec.compParams(26);
            reec.pfflag = reec.compParams(27);
            reec.vflag = reec.compParams(28);
            reec.qflag = reec.compParams(29);
            reec.pqflag = reec.compParams(30);
            reec.stWrite=fopen(reec.reecOutputLocation,'w');
            if reec.stWrite<0
                sprintf('error opening file')
            else
                fprintf(reec.stWrite,'%8s %8s %8s %8s %8s %8s \r\n',['State 1' 'State 2' 'State 3' 'State 4' 'State 5' 'State 6']);
            end  
        end
        %% Writing exciter states    
        function reec=WriteData(reec)
                fprintf(reec.stWrite,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \r\n',[reec.reecStates(1) reec.reecStates(2) reec.reecStates(3) reec.reecStates(4) reec.reecStates(5) reec.reecStates(6)]);  
        end
        %%
        function [reec,qext,pref]=Init(reec,ipcmd,iqcmd,vterm,qgen,pe)
            ip=ipcmd;
            pref=vterm*ip;
            reec.vref0=vterm;
            vtfilt=vterm;
            
            if reec.qflag==0
                qout1=iqcmd;
                q51=0;
                q31=0;
                qin1=qout1;
                qin=vterm*qin1;
            elseif reec.qflag==1
                qout1=0;
                q51=iqcmd;
                q34=vterm;
                if reec.vflag==1
                    q3=q34;
                    q31=q3;
                    qin=qgen;
                elseif reec.vflag==0
                    qin=q34;
                    q31=0;
                end
            end
            
            if reec.pfflag==0
                qext=qin;
                reec.pfaref=0;
                p111=0;
            elseif reec.pfflag==1
                reec.pfaref=atan(qin/pe);
                p111=pe;
            end
            reec.reecStates=[vtfilt;p111;q31;q51;qout1;pref]; % 6 X 1 vector of initial conditions
        end
        %%
        function [reec,iqcmd,ipcmd]=newStates(reec,tspan,pe,vt,qext,qgen,pref,ipcmd0,iqcmd0)
            options=odeset('RelTol',1e-9);
          
             if reec.pqflag==0
                 iqmax=reec.imax;
                 iqmin=-iqmax;
                 ipmax=sqrt(reec.imax^2-iqcmd0^2);
                 ipmin=-ipmax; % made this change to allow negative generation for storage
             elseif reec.pqflag==1
                 iqmax=sqrt(reec.imax^2-ipcmd0^2);
                 iqmin=-iqmax;
                 ipmax=reec.imax;
                 ipmin=-ipmax; % made this change to allow negative generation for storage
             end
            
            [T,X]=ode23tb(@(t,x) reec_dyn(t,x,pe,vt,qext,qgen,pref,reec.vdip,reec.vup,reec.trv,reec.pfaref,reec.tp,reec.qmax,reec.qmin,reec.vmax,reec.vmin,reec.kqp,reec.kqi,reec.kvp,reec.kvi,reec.tiq,reec.dpmax,reec.dpmin,reec.pmax,reec.pmin,reec.tpord,reec.pfflag,reec.vflag,reec.qflag),tspan,reec.reecStates,options);
            
            function dy=reec_dyn(t1,y,pe,vt,qext,qgen,pref,vdip,vup,trv,pfaref,tp,qmax,qmin,vmax,vmin,kqp,kqi,kvp,kvi,tiq,dpmax,dpmin,pmax,pmin,tpord,pfflag,vflag,qflag)
                dy=zeros(6,1); 
                            % Voltage Transducer Dynamics
                if trv==0
                    dy(1)=0;
                    vtfilt=vt;
                else
                    dy(1)=(vt-y(1))/trv;
                    vtfilt=y(1);
                end
                                
                if pfflag==1
                    tt=tan(pfaref);
                    dy(2)=(pe-y(2))/tp;
                    qin=y(2)*tt;
                else
                    dy(2)=0;
                    qin=qext;
                end
                
                if qflag==0
                    vt1=vtfilt;
                    if vt1<=0.01
                        vt1=0.01;
                    end
                    qin1=qin/vt1;
                    dy(5)=(qin1-y(5))/tiq;
                    dy(3)=0;
                    dy(4)=0;
                elseif qflag==1
                    q1=qin;
                    if q1>=qmax
                        q1=qmax;
                    elseif q1<=qmin
                        q1=qmin;
                    end
                    
                    q2=q1-qgen;
                    q3=kqp*q2+y(3);
                    dy(3)=kqi*q2;
                    if ((q3<=vmin && dy(3)<0))
                        dy(3)=0;
                    elseif (q3>=vmax && dy(3)>0)
                        dy(3)=0;
                        if dy(4)>0
                            dy(4)=0;
                        end
                    end
                    if vflag==1
                        q34=q3;
                    elseif vflag==0
                        q34=qin;
                    end
                    if q34>=vmax
                        q34=vmax;
                    elseif q34<=vmin
                        q34=vmin;
                    end
                    q4=q34-vtfilt;
                    q5=kvp*q4+y(4);
                    dy(4)=kvi*q4;
                    dy(5)=0;
                    if q5<=iqmin && dy(4)<0
                        dy(4)=0;
                        if dy(3)<0
                            dy(3)=0;
                        end
                    elseif q5>=iqmax && dy(4)>0
                        dy(4)=0;
                    end
                end
                if tpord==0
                    dy(6)=0;
                else
                    dy(6)=(pref-y(6))/tpord;
                    if y(6)>=pmax && dy(6)>0
                        dy(6)=0;
                    elseif y(6)<=pmin && dy(6)<0
                        dy(6)=0;
                    elseif y(6)<pmax && y(6)>pmin && dy(6)>=dpmax
                        dy(6)=dpmax;
                    elseif y(6)<pmax && y(6)>pmin && dy(6)<=dpmin
                        dy(6)=dpmin;
                    end
                end
                if vt<vdip ||vt>vup
                    dy(3:6)=0;
                end
                if norm(dy)>1e-5
                    xxxx=55;
                end
            end
            reec.reecStates=X(end,:);
            
            %%
                
                if reec.trv==0
                    y1=-vt+reec.vref0;
                    vtfilt_1=vt;
                    reec.reecStates(1)=vt;
                else
                    y1=-reec.reecStates(1)+reec.vref0;
                    vtfilt_1=reec.reecStates(1);
                end
                if (y1<=0 && y1>=reec.dbd1) || (y1>0 && y1<=reec.dbd2)
                    y2=0;
                else
                    y2=y1;
                end
                y3=reec.kqv*y2;
                y4=y3;
                if y3>=reec.iqhl
                    y4=reec.iqhl;
                elseif y3<=reec.iqll
                    y4=reec.iqll;
                end
                iqinj=y4;
                                
                if reec.pfflag==1
                    tt_1=tan(reec.pfaref);
                    qin_1=reec.reecStates(2)*tt_1;
                else
                    qin_1=qext;
                end
                
                if reec.qflag==0
                    q6=reec.reecStates(5);
                elseif reec.qflag==1
                    q1_1=qin_1;
                    if q1_1>=reec.qmax
                        q1_1=reec.qmax;
                    elseif q1_1<=reec.qmin
                        q1_1=reec.qmin;
                    end
                    
                    q2_1=q1_1-qgen;
                    q3_1=reec.kqp*q2_1+reec.reecStates(3);
                    if reec.vflag==1
                        q34_1=q3_1;
                    elseif reec.vflag==0
                        q34_1=qin_1;
                    end
                    if q34_1>=reec.vmax
                        q34_1=reec.vmax;
                    elseif q34_1<=reec.vmin
                        q34_1=reec.vmin;
                    end
                    q4_1=q34_1-vtfilt_1;
                    q5_1=reec.kvp*q4_1+reec.reecStates(4);
                    q6=q5_1;
                end
                iqcmd=iqinj+q6;
                if iqcmd>=iqmax
                    iqcmd=iqmax;
                elseif iqcmd<=iqmin
                    iqcmd=iqmin;
                end
               
                vt1_1=vtfilt_1;
                if vt1_1<=0.01
                    vt1_1=0.01;
                end
                if (reec.tpord==0)
                    ipcmd=pref/vt1_1;
                else
                    ipcmd=reec.reecStates(6)/vt1_1;
                end
                if ipcmd>=ipmax
                    ipcmd=ipmax;
                elseif ipcmd<=ipmin
                    ipcmd=ipmin;
                end
            
        end
        
    end
    
end

