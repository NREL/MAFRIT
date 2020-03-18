classdef reec_a < Components
    % REEC_A represents the electrical control module for Type 3 and 4 Wind Turbines and Central PV Stations.
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
        iqfrz;
        thld;
        thld2;
        pfaref; %this parameter is not input, but is assigned during initialization
        tp;
        qmax;
        qmin;
        vmax;
        vmin;
        kqp;
        kqi;
        kvp;
        kvi;
        vref1;
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
        pflag;
        pqflag;
        vq1;
        iq1;
        vq2;
        iq2;
        vq3;
        iq3;
        vq4;
        iq4;
        vp1;
        ip1;
        vp2;
        ip2;
        vp3;
        ip3;
        vp4;
        ip4;
        voltage_dip;
        iqinj;
        reset;
        ipcmd_old;
        reecStates;
        reecOutputLocation;
        stWrite;
    end
    
    methods
        function reec=reec_a()
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
            reec.iqfrz = reec.compParams(11);
            reec.thld = reec.compParams(12);
            reec.thld2 = reec.compParams(13);
            reec.tp = reec.compParams(14);
            reec.qmax = reec.compParams(15);
            reec.qmin = reec.compParams(16);
            reec.vmax = reec.compParams(17);
            reec.vmin = reec.compParams(18);
            reec.kqp = reec.compParams(19);
            reec.kqi = reec.compParams(20);
            reec.kvp = reec.compParams(21);
            reec.kvi = reec.compParams(22);
            reec.vref1 = reec.compParams(23);
            reec.tiq = reec.compParams(24);
            reec.dpmax = reec.compParams(25);
            reec.dpmin = reec.compParams(26);
            reec.pmax = reec.compParams(27);
            reec.pmin = reec.compParams(28);
            reec.imax = reec.compParams(29);
            reec.tpord = reec.compParams(30);
            reec.pfflag = reec.compParams(31);
            reec.vflag = reec.compParams(32);
            reec.qflag = reec.compParams(33);
            reec.pflag = reec.compParams(34);
            reec.pqflag = reec.compParams(35);
            reec.vq1 = reec.compParams(36);
            reec.iq1 = reec.compParams(37);
            reec.vq2 = reec.compParams(38);
            reec.iq2 = reec.compParams(39);
            reec.vq3 = reec.compParams(40);
            reec.iq3 = reec.compParams(41);
            reec.vq4 = reec.compParams(42);
            reec.iq4 = reec.compParams(43);
            reec.vp1 = reec.compParams(44);
            reec.ip1 = reec.compParams(45);
            reec.vp2 = reec.compParams(46);
            reec.ip2 = reec.compParams(47);
            reec.vp3 = reec.compParams(48);
            reec.ip3 = reec.compParams(49);
            reec.vp4 = reec.compParams(50);
            reec.ip4 = reec.compParams(51);
            
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
        function [reec,pord,qext,pref]=Init(reec,ipcmd,iqcmd,vterm,wg,qgen,pe)
            ip=ipcmd;
            pord=vterm*ip;
            reec.vref0=vterm;
            vtfilt=vterm;
            p2=pord;
            if reec.pflag==1
                p1=p2/wg;    % Since Type 4 wind turbine models can only include wtgt model, w0 must be user defined unlike type 3 turbines where w0=wref as calculated from the piecewise linear function f(pe) in wtgtq model. The appropriate value of w0 should be calculated in the main program itself.        
            else
                p1=p2;
            end
            pref=p1;
            
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
                    qin=q34-reec.vref1;
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
            
            reec.voltage_dip=0;
            reec.reset=0;
            reec.iqinj=0;
            reec.reecStates=[vtfilt;p111;q31;q51;qout1;pord]; % 6 X 1 vector of initial conditions
        end
        %%
        function [reec,iqcmd,ipcmd,pord]=newStates(reec,tspan,pe,wg,vt,qext,qgen,pref,ipcmd0,iqcmd0)
            options=odeset('RelTol',1e-9);
            deltat=tspan(2)-tspan(1);
            persistent cnt
            if isempty(cnt)
                cnt=0.00001;
            end
            
            % voltage dip logic
            if vt<reec.vdip || vt>reec.vup
                reec.voltage_dip=1;  
            else
                vdip_prev=reec.voltage_dip;
                reec.voltage_dip=0;
                if vdip_prev>reec.voltage_dip
                    reec.reset=1;
                end   
            end
            
            if reec.reecStates(1)>=0 && reec.reecStates(1)<reec.vq1
                iq=reec.iq1;
            elseif reec.reecStates(1)>=reec.vq1 && reec.reecStates(1)<reec.vq2
                iq=reec.iq1+(reec.iq2-reec.iq1)*(reec.reecStates(1)-reec.vq1)/(reec.vq2-reec.vq1);
            elseif reec.reecStates(1)>=reec.vq2 && reec.reecStates(1)<reec.vq3
                iq=reec.iq2+(reec.iq3-reec.iq2)*(reec.reecStates(1)-reec.vq2)/(reec.vq3-reec.vq2);
            elseif reec.reecStates(1)>=reec.vq3 && reec.reecStates(1)<reec.vq4
                iq=reec.iq3+(reec.iq4-reec.iq3)*(reec.reecStates(1)-reec.vq3)/(reec.vq4-reec.vq3);
            else
                iq=reec.iq4;
            end
            
            if reec.reecStates(1)>=0 && reec.reecStates(1)<reec.vp1
                ip=reec.ip1;
            elseif reec.reecStates(1)>=reec.vp1 && reec.reecStates(1)<reec.vp2
                ip=reec.ip1+(reec.ip2-reec.ip1)*(reec.reecStates(1)-reec.vp1)/(reec.vp2-reec.vp1);
            elseif reec.reecStates(1)>=reec.vp2 && reec.reecStates(1)<reec.vp3
                ip=reec.ip2+(reec.ip3-reec.ip2)*(reec.reecStates(1)-reec.vp2)/(reec.vp3-reec.vp2);
            elseif reec.reecStates(1)>=reec.vp3 && reec.reecStates(1)<reec.vp4
                ip=reec.ip3+(reec.ip4-reec.ip3)*(reec.reecStates(1)-reec.vp3)/(reec.vp4-reec.vp3);
            else
                ip=reec.ip4;
            end
             if reec.pqflag==0
                 iqmax=min(iq,reec.imax);
                 iqmin=-iqmax;
                 ipmax=min(ip,sqrt(reec.imax^2-iqcmd0^2));
                 ipmin=0;
             elseif reec.pqflag==1
                 iqmax=min(iq,sqrt(reec.imax^2-ipcmd0^2));
                 iqmin=-iqmax;
                 ipmax=min(ip,reec.imax);
                 ipmin=0;
             end
            
            [T,X]=ode23t(@(t,x) reec_dyn(t,x,pe,wg,vt,qext,qgen,pref,reec.vdip,reec.vup,reec.trv,reec.pfaref,reec.tp,reec.qmax,reec.qmin,reec.vmax,reec.vmin,reec.kqp,reec.kqi,reec.kvp,reec.kvi,reec.vref1,reec.tiq,reec.dpmax,reec.dpmin,reec.pmax,reec.pmin,reec.tpord,reec.pfflag,reec.vflag,reec.pflag,reec.qflag),tspan,reec.reecStates,options);
            
            function dy=reec_dyn(t1,y,pe,wg,vt,qext,qgen,pref,vdip,vup,trv,pfaref,tp,qmax,qmin,vmax,vmin,kqp,kqi,kvp,kvi,vref1,tiq,dpmax,dpmin,pmax,pmin,tpord,pfflag,vflag,pflag,qflag)
                dy=zeros(6,1); 
                        % Voltage Transducer Dynamics
                if trv==0
                    dy(1)=0;
                    vtfilt=vt;
                else
                    dy(1)=(vt-y(1))/trv;
                    vtfilt=y(1);
                end
%                         % Voltage_Dip Logic
%                 if vt<vdip ||vt>vup
%                     reec.voltage_dip=1;
%                 else
%                     reec.voltage_dip=0;
%                 end
                        % Reactive Power Control        
                if pfflag==1
                    tt=tan(pfaref);
                    dy(2)=(pe-y(2))/tp;
                    qin=y(2)*tt;
                else
                    dy(2)=0;
                    qin=qext;
                end
                        % Reactive Control and Local Voltage Response
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
                                q34=qin+vref1;
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
                        % Active Power Control
                        if reec.tpord~=0
                            p1=pref;
                            if pflag==1
                                p2=p1*wg;
                            elseif pflag==0
                                p2=p1;
                            end
                            dy(6)=(p2-y(6))/tpord;
                            if y(6)>=pmax && dy(6)>0
                                dy(6)=0;
                            elseif y(6)<=pmin && dy(6)<0
                                dy(6)=0;
                            end
                        end
                if reec.voltage_dip==1
                    dy(3:6)=0;
                end
            end
            reec.reecStates=X(end,:);  
            
            %% 
                if reec.trv==0
                    y1_1=-vt+reec.vref0;
                    vtfilt_1=vt;
                    reec.reecStates(1)=vt;
                else
                    y1_1=-reec.reecStates(1)+reec.vref0;
                    vtfilt_1=reec.reecStates(1);
                end 
                if (y1_1<=0 && y1_1>=reec.dbd1) || (y1_1>0 && y1_1<=reec.dbd2)
                    y2_2=0;
                else
                    y2_2=y1_1;
                end
                y3_3=reec.kqv*y2_2;
                y4_4=y3_3;
                if y3_3>=reec.iqhl
                    y4_4=reec.iqhl;
                elseif y3_3<=reec.iqll
                    y4_4=reec.iqll;
                end
                
                if (reec.voltage_dip==1 && reec.reset==0)
                    reec.iqinj=y4_4;
                elseif (reec.reset==1 && reec.thld==0) || (reec.voltage_dip==0 && reec.thld~=0 && cnt>abs(reec.thld))
                    reec.iqinj=0;
                    cnt=0.00001;
                    reec.reset=0;
                elseif (reec.reset==1 && reec.thld>0 && cnt<abs(reec.thld))
                    reec.iqinj=y4_4;
                    cnt=cnt+deltat;
                elseif (reec.reset==1 && reec.thld<0 && cnt<abs(reec.thld))
                    reec.iqinj=reec.iqfrz;
                    cnt=cnt+deltat;
                end
                
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
                    
                    q2_2=q1_1-qgen;
                    q3_3=reec.kqp*q2_2+reec.reecStates(3);
                    if reec.vflag==1
                        q34_1=q3_3;
                    elseif reec.vflag==0
                        q34_1=qin_1+reec.vref1;
                    end
                    if q34_1>=reec.vmax
                        q34_1=reec.vmax;
                    elseif q34_1<=reec.vmin
                        q34_1=reec.vmin;
                    end
                    q4_4=q34_1-vtfilt_1;
                    q5_1=reec.kvp*q4_4+reec.reecStates(4);
                    q6=q5_1;
                end
                iqcmd=reec.iqinj+q6;
                if iqcmd>=iqmax
                    iqcmd=iqmax;
                elseif iqcmd<=iqmin
                    iqcmd=iqmin;
                end
                
                pord=reec.reecStates(6);
                vt1_1=vtfilt_1;
                if vt1_1<=0.01
                    vt1_1=0.01;
                end
                
                if (reec.reset==0 && reec.tpord~=0)
                    ipcmd=reec.reecStates(6)/vt1_1;
                    reec.ipcmd_old=ipcmd;
                elseif (reec.reset==0 && reec.tpord==0)
                    p1=pref;
                    if reec.pflag==1
                        p2=p1*wg;
                    elseif reec.pflag==0
                        p2=p1;
                    end
                    ipcmd=p2/vt1_1;
                    reec.ipcmd_old=ipcmd;
                elseif (reec.reset==1)
                    ipcmd=reec.ipcmd_old;
                end
                if ipcmd>=ipmax
                    ipcmd=ipmax;
                elseif ipcmd<=ipmin
                    ipcmd=ipmin;
                end
            
            
        end
        
    end
    
end

