classdef ggov1 < Components
    % GGOV1 Initializes all governors of type GGOV1 and updates governor
    % states during simulation. It is based on the block diagram given here:
    % https://www.nerc.com/comm/PC/NERCModelingNotifications/Gas_Turbine_Governor_Modeling.pdf
    
    properties
        sbturb;
        r;
        rselect;
        tpelec;
        maxerr;
        minerr;
        kpgov;
        kigov;
        kdgov;
        tdgov;
        vmax;
        vmin;
        tact;
        kturb;
        wfnl;
        tb;
        tc;
        flag;
        teng;
        tfload;
        kpload;
        kiload;
        ldref;
        dm;
        ropen;
        rclose;
        kimw;
        pmwset;
        aset;
        ka;
        ta;
        db;
        tsa;
        tsb;
        rup;
        rdown;
        fsr1;               % fsr1 and cfe1 are not governor parameters, but inputs into the dynamics blocks of ggov1. These had to defined as governor properties to break the algebraic loop in the ggov1 model. 
        cfe1;
        govStates;
        govOutputLocation;
        stWrite;
    end
    
    methods
        function gov=ggov1()
        end
        
        function gov=AssignParams(gov)
            gov.sbturb=gov.compParams(1);
            gov.r=gov.compParams(2);
            gov.rselect=gov.compParams(3);
            gov.tpelec=gov.compParams(4);
            gov.maxerr=gov.compParams(5);
            gov.minerr=gov.compParams(6);
            gov.kpgov=gov.compParams(7);
            gov.kigov=gov.compParams(8);
            gov.kdgov=gov.compParams(9);
            gov.tdgov=gov.compParams(10);
            gov.vmax=gov.compParams(11);
            gov.vmin=gov.compParams(12);
            gov.tact=gov.compParams(13);
            gov.kturb=gov.compParams(14);
            gov.wfnl=gov.compParams(15);
            gov.tb=gov.compParams(16);
            gov.tc=gov.compParams(17);
            gov.flag=gov.compParams(18);
            gov.teng=gov.compParams(19);
            gov.tfload=gov.compParams(20);
            gov.kpload=gov.compParams(21);
            gov.kiload=gov.compParams(22);
            gov.ldref=gov.compParams(23);
            gov.dm=gov.compParams(24);
            gov.ropen=gov.compParams(25);
            gov.rclose=gov.compParams(26);
            gov.kimw=gov.compParams(27);
            gov.pmwset=gov.compParams(28);
            gov.aset=gov.compParams(29);
            gov.ka=gov.compParams(30);
            gov.ta=gov.compParams(31);
            gov.db=gov.compParams(32);
            gov.tsa=gov.compParams(33);
            gov.tsb=gov.compParams(34);
            gov.rup=gov.compParams(35);
            gov.rdown=gov.compParams(36);
            gov.stWrite=fopen(gov.govOutputLocation,'w');
            if gov.stWrite<0
                sprintf('error opening file')
            else
                fprintf(gov.stWrite,'%8s %8s %8s %8s %8s %8s %8s %8s %8s %8s \r\n',['State 1' 'State 2' 'State 3' 'State 4' 'State 5' 'State 6' 'State 7' 'State 8' 'State 9' 'State 10']);
            end 
        end
%%        
        function [gov,pref]=Init(gov,sbgen,pelec)
            pelect=pelec*sbgen/gov.sbturb; 
            
            if gov.dm<=0
                y11=pelect;
            else
                y11=pelect+gov.dm;
            end
            
            u11=y11/gov.kturb;
            cfe=u11-gov.wfnl;
            y10=cfe;                                % State 9
            fsr=y10;
            y111=(1-gov.tc/gov.tb)*u11*gov.kturb;   % State 10
            if gov.tpelec==0
                y1=0;                               % State 1
            else
                y1=pelect;                          % State 1 
            end
            gov.pmwset=y1;
            
            switch gov.rselect
                case 1
                    u3=1+y1*gov.r;
                case -2
                    u3=1+gov.r*fsr;
                case -1
                    u3=1+gov.r*y10;
                case 0
                    u3=1;
                otherwise
                    sprintf('%s','non-feasible value of rselect; check ggov1 parameters')
                    pause
            end
            
            y431=0;                                 % State 4
            u51=1;                                  % State 5
            if gov.dm<0
                u7=gov.dm*y10;
            else
                u7=y10;
            end
            
            y6=u7;
            y7=y6;                                  % State 7
            u6=u7;
            if gov.tsb==0
                y61=0;                              % State 6
            else
                y61=u7*(1-gov.tsa/gov.tsb);         % State 6
            end
            
            gov.ldref=(u6-gov.wfnl)*gov.kturb;
            
            y4=y10; % y4 is fsrn
            y42=y4;                                 % State 3
            y8=y4;  % y8 is fsrt
            y81=y8;                                 % State 8
            
            if gov.kimw==0
                y2=0;
                pref=u3;  %pref=u3 if kimw=0; pref cannot be chosen independently when kimw=0;
            else
                pref=1; % if kimw~=0 then pref can be chosen independently because at steady state pref+y2=u3 and there is no relation that can help us obtain unique values of pref and y2.
                y2=u3-pref;                          % State 2
            end
            gov.govStates=[y1;y2;y42;y431;u51;y61;y7;y81;y10;y111];
            gov.fsr1=y10;
            gov.cfe1=y10;
        end
        %% Writing governor states    
        function gov=WriteData(gov)
                fprintf(gov.stWrite,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\r\n',[gov.govStates(1) gov.govStates(2) gov.govStates(3) gov.govStates(4) gov.govStates(5) gov.govStates(6) gov.govStates(7) gov.govStates(8) gov.govStates(9) gov.govStates(10)]);    
        end
        %%
        function [gov,pmech]=newStates(gov,tspan,pref,pelec,wr)
            options=odeset('RelTol',1e-9);
            fsr=gov.fsr1;  % fsr and cfe span the ode function and the function containing the derivatives. Therefore, any change in these variables during integration is reflected simultaneously outside the ode function as well. 
            cfe=gov.cfe1;
            
            [T,X]=ode23t(@(t,x) gov_dyn(t,x,pref,pelec,wr,gov.r,gov.rselect,gov.tpelec,gov.maxerr,gov.minerr,gov.kpgov,gov.kigov,gov.kdgov,gov.tdgov,gov.vmax,gov.vmin,gov.tact,gov.kturb,gov.wfnl,gov.tb,gov.tc,gov.flag,gov.teng,gov.tfload,gov.kpload,gov.kiload,gov.ldref,gov.dm,gov.ropen,gov.rclose,gov.kimw,gov.pmwset,gov.aset,gov.ka,gov.ta,gov.db,gov.tsa,gov.tsb,gov.rup,gov.rdown),tspan,gov.govStates,options);
            
            function dy=gov_dyn(t1,y,pref,pelec,wr,r,rselect,tpelec,maxerr,minerr,kpgov,kigov,kdgov,tdgov,vmax,vmin,tact,kturb,wfnl,tb,tc,flag,teng,tfload,kpload,kiload,ldref,dm,ropen,rclose,kimw,pmwset,aset,ka,ta,db,tsa,tsb,rup,rdown)
                
                dy=zeros(10,1);
                
                % Supervisory Block Dynamics
                
                    if tpelec==0
                        dy(1)=0;
                    else
                        dy(1)=(pelec-y(1))/tpelec;
                    end

                    if kimw==0
                        dy(2)=0;
                    else
                        dy(2)=kimw*(pmwset-y(1));
                    end

                    if y(2)>=1.1*r || y(2)<=-1.1*r
                        dy(2)=0;
                    end             

                    u3=pref+y(2);
                    switch rselect
                        case 1
                            x=y(1);
                        case -2
                            x=1+r*fsr;  
                        case -1
                            x=y(9);
                        case 0
                            x=0;
                    end

                    y3=u3-x*r-wr;


                    if y3>=db || y3<=-db
                        u4=y3;
                    else
                        u4=0;
                    end

                    if u4>=maxerr
                        u4=maxerr;
                    elseif u4<=minerr
                        u4=minerr;
                    end
                
                % PID Governor Dynamics

                    y41=kpgov*u4;
                    dy(3)=kigov*u4;
                    if kdgov==0
                        dy(4)=0;
                        y43=0;
                    else
                        dy(4)=((kdgov/tdgov)*u4-y(4))/tdgov;
                        y43=u4*(kdgov/tdgov)-y(4);
                    end
                    y4=y41+y(3)+y43;
                    fsrn=y4;
                
                % Speed and Acceleration Effects
                    if ka==0
                        y5=9999;
                        dy(5)=0;
                    else
                        dy(5)=(wr-y(5))/ta;
                        u5=wr/ta-y(5);
                        y5=(aset-u5)*ka+fsr; 
                    end
                    fsra=y5;
               
                 % Load Reference Controller
                    if dm<0
                        u7=wr*dm*cfe;   
                    else
                        u7=cfe;
                    end
                    
                    if tsb==0
                        y6=u7;
                        dy(6)=0;
                    else
                        dy(6)=((1-tsa/tsb)*u7-y(6))/tsb;
                        y6=(tsa/tsb)*u7+y(6);
                    end
                    
                   dy(7)=(y6-y(7))/tfload;
                   
                   u6=ldref/kturb+wfnl;
                   u8=u6-y(7);
                   dy(8)=kiload*u8;
                   y8=kpload*u8+y(8);
                
                   if y8>=rup
                       fsrt=rup;
                   elseif y8<=rdown
                       fsrt=rdown;
                   else
                       fsrt=y8;
                   end
                   
                % Low Value Select Logic   
                    y9=min([fsrn,fsra,fsrt]);
                    
                    if y9>=vmax
                        fsr=vmax;
                    elseif y9<=vmin
                        fsr=vmin;
                    else
                        fsr=y9;
                    end
                    
                %  Gate Valve Actuator Dynamics
                
                    dy(9)=(fsr-y(9))/tact;
                    if (y(9)>=ropen && dy(9)>0) || (y(9)<=rclose && dy(9)<0)
                        dy(9)=0;
                    end
                
                %  Turbine Dynamics
                
                    if flag==1
                        cfe=y(9);
                    else
                        cfe=y(9)*wr;
                    end
                    
                    u11=cfe+wfnl;
                    
                    dy(10)=((1-tc/tb)*u11*kturb-y(10))/tb;
                                     
            end
            gov.govStates=X(end,:);  
            gov.fsr1=fsr;
            gov.cfe1=cfe;
            if gov.flag==1
                cfe11= gov.govStates(9);
            else
                cfe11= gov.govStates(9)*wr;
            end
            u111=cfe11+gov.wfnl;
            y11=(gov.tc/gov.tb)*u111*gov.kturb+gov.govStates(10);
               
            if gov.dm>0
                pmech=y11-gov.dm*wr;
            else
                pmech=y11;          % pmech is not a direct output of a governor state, but a function of multiple states. 
                                    % Therefore, it is calculated based on the last state determined by MATLAB
            end

        end
        
    end
    
end

